function [] =  EstimatePL_realdata(loopNum,dateStart,dataName,numSamples,paramToEstimate,addRandomness,increasebound,lowersteps,uppersteps,refine,smooth)
% Function to run profile likelihood for real data, starting at n=numSamples
% fixed steps, where each loopNum=1:numSamples optimizes one step.

%% Setup
s=rng('shuffle'); %do not note remove this line
d = datestr(now,'YYmmDD-HHMMSSFFF'); %creating different seed for each start
s=rng(s.Seed+loopNum*sum(d)); %do not note remove this line

basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');

%add all project folders to the matlab path
addpath(genpath(fullfile(basefolder,'Optimization')))
addpath(genpath(fullfile(basefolder,'Data')))
addpath(genpath(fullfile(basefolder,'Modelfiles')))
addpath(genpath(fullfile(basefolder,'Requirements')))
addpath(genpath(fullfile(basefolder,'Simulation')))
addpath(genpath(fullfile(basefolder,'Tools')))

% setup AMICI and MEIGO toolboxes
run(fullfile(basefolder, 'Requirements', 'AMICI-0.10.11_SS_eventFix', 'matlab', 'installAMICI.m'))
run(fullfile(basefolder, 'Requirements', 'MEIGO', 'install_MEIGO.m'))


if isempty(paramToEstimate)
    paramNamesToEstimate = {'m2_LV' 'Caa' 'Rao'};
else
    paramNamesToEstimate = {paramToEstimate};
end

if nargin < 11
    smooth =0;
    refine = 0;
end

if nargin < 9
    lowersteps = 0;
    uppersteps = 0;
end

if nargin < 7
    increasebound = 0;
end

if nargin < 6
    addRandomness = 0;
end

addRandomnessSmall = 0;


%% Load data
load(fullfile(basefolder,'Data',[dataName, '.mat']),'estimationData')
data = estimationData;
disp(dataName)

%% Load parameters and data
resultsfolder = fullfile(basefolder,'Parameters'); 
datain.(dataName) = data;
[~,data,~,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues,~,~,~,~,~] = setup_simulations({dataName},datain,resultsfolder,1);

% settings for the cost function:
doPlot = 0;
dispErrors = 0;

%% number of optimizations for each parameter value
if smooth && refine %add extra steps close to limits of the confidence interval
    nRepeats = 30;
elseif addRandomness
    nRepeats = 40+10;%7   %5 extra in the end to improve the best found (10 last reps are not using random startguess)
else
    nRepeats = 40;%7
end

%% ESS options
% MEIGO OPTIONS I (COMMON TO ALL SOLVERS):
opts.ndiverse     = 'auto'; %100; %500; %5; %
if strcmp(data.meta.simulationFunction,'simulate_avatar_correction')
    opts.maxtime      = 2*600; % MAX-Time of optmization, i.e how long the optimization will last %Maximum CPU time in seconds (Default 60)
else
    opts.maxtime      = 600; % MAX-Time of optmization, i.e how long the optimization will last %Maximum CPU time in seconds (Default 60)
end
opts.maxeval      = 1e8; % max number of evals, i.e cost function calls (Default 1000)
opts.log_var      = [];  %skip this

opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; % local solver, fmincon works in my experience best
opts.local.finish = opts.local.solver; %uses the local solver to check the best p-vector
opts.local.bestx = 0; % read the documentation, think it's best to leave at zero for now.
problem.f   = 'costFunction_chi2'; % %name of cost-function
opts.iterprint = 0;

% MEIGO OPTIONS II (FOR ESS AND MULTISTART):
opts.local.iterprint = 0; % prints what going on during optimization

% MEIGO OPTIONS III (FOR ESS ONLY):
opts.dim_refset   = 'auto'; % leave to auto for now

% OPTIONS AUTOMATICALLY SET AS A RESULT OF PREVIOUS OPTIONS:
if(strcmp(opts.local.solver,'fmincon'))
    opts.local.use_gradient_for_finish = 1; %DW: provide gradient to fmincon
else
    opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
end
opts.local.check_gradient_for_finish = 0; %DW: gradient checker


%% FIND STARTGUESS AT OPTIMUM
folders = dir(fullfile(resultsfolder,'ESS',['E_' dataName '*']));

if isempty(folders)
    disp('ERROR: no startguess found')
else
    [~,latestfolderInd] = max([folders.datenum]);
    foldernames = {folders.name};
    folderpaths = {folders.folder};
    loadfolderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    [bestParams,optcost,~,loadedconstants,~,~,~,optcostChi2] = findBestParams(loadfolderName,0,length(paramNames));
    fprintf('Loaded paramvalues with cost %0.2f\n',optcostChi2)
end

%% Set bounds (separate for each subject, since paramdata is individual)
[lbOrig,ubOrig] = loadParamBounds(ind,data.parameters);
lbOrig = lbOrig(data.meta.paramsToOptimize);
ubOrig = ubOrig(data.meta.paramsToOptimize);

allparams = origParamvalues;

%% find folder to save results in for this data
PLfolder = fullfile(basefolder,'Parameters','PL',['E_' dataName]);

%% Run optimization
format long
format compact
warning('off','all') % AMICI prints error-messages for all simulations (if they fail) so this will fix that annoying orange text form appearing

for n = 1:length(paramNamesToEstimate)
    paramname = paramNamesToEstimate{n};
    indFixedParam = strcmp(paramname,paramNames);
    fprintf('Starting paramestimation of %s for %s \n',paramname,dataName)

    %% Set test ranges
    %load bounds found from sampling
    load(fullfile(basefolder,'Optimization','optbounds.mat'),'paramnames','lb','ub')
    indparam = strcmp(paramname,paramnames);
    minval = lb(indparam) - lb(indparam)*0.5;% 0.2 was too small, should increase.
    maxval = ub(indparam) + ub(indparam)*1;% 0.2 was too small, should increase.
    optval = bestParams(indFixedParam);

    if optval <= minval
        minval = max(0,optval  - optval*0.4);
    elseif optval >= maxval
        maxval = optval  + optval*0.4;
    end
    testrange = linspace(minval,maxval,numSamples);

    %% Load previously optmized range if needed
    if increasebound || smooth || refine
        % load all PLs for this experiment
        loadParameters=1;
        [~,allPLs,~,~,bestParams,bestPLCost] = findPLparams(dataName,resultsfolder,length(paramNames),loadParameters,{paramname},paramNames);
        optval = bestParams(indFixedParam);
        if length(bestParams) > 28
            bestParams = bestParams(1:28);
        end

        %sort
        [~,sind]=sort(allPLs.(paramname)(:,1));

        %check confidence interval range
        costs = allPLs.(paramname)(sind,2);
        pvals = allPLs.(paramname)(sind,1);
        bestcost = min(costs);
        limit = bestcost+chi2inv(0.95,1);
        withinrange = pvals(costs<=limit);
        lbparam = min(withinrange);
        ubparam = max(withinrange);
    end

    %% Add more steps closer to the chi2 limit
    if refine && smooth
        disp('reine+smooth = add 4 steps closer to parameter chi2 limits')
        highervals = pvals(pvals > ubparam)';
        lowervals = pvals(pvals < lbparam)';
        lowervals = sort(lowervals, 'descend');
        highervals = sort(highervals, 'ascend');
        if isempty(lowervals)
            lowervals = lbparam*0.6;
        end
        if isempty(highervals)
            highervals = ubparam*1.4;
        end
        diffh = min(highervals)-ubparam;
        diffl = lbparam-max(lowervals);
        %add steps between the limit values and the next value that is not
        %below the limit
        testrange = [ubparam+4*diffh/5,ubparam+3*diffh/5,ubparam+2*diffh/5,ubparam+diffh/5,lbparam-4*diffl/5,lbparam-3*diffl/5,lbparam-2*diffl/5,lbparam-diffl/5];
        testrange = repmat(testrange,1,ceil(numSamples/length(testrange)));% --> range is at least numSamples long (ie if is 32 --> 4 repeats)
        testrange = testrange(1:numSamples);%--> testrange is exactly numSamples long

        %% Refine: add more steps to refine the parameter limits
    elseif refine
        disp('Refine bounds')
        highervals = pvals(pvals > ubparam)';
        lowervals = pvals(pvals < lbparam)';
        lowervals = sort(lowervals, 'descend');
        highervals = sort(highervals, 'ascend');
        if ubparam == max(pvals)
            %increase range
            testrangeU = linspace(ubparam,ubparam*5,ceil(numSamples/2));
        else
            testrangeU = highervals(1:min(ceil(numSamples/2),length(highervals)));
            if length(testrangeU) < ceil(numSamples/2)
                %add extra steps to have exactly numSamples/2 numbers in
                %the range
                l = length(testrangeU);
                testrangeU = [testrangeU,linspace(ubparam,ubparam*5,ceil(numSamples/2)-l)];
            end
        end

        if lbparam == min(pvals)
            %increase range
            testrangeL = linspace(lbparam/5,lbparam,floor(numSamples/2));
        else
            testrangeL = lowervals(1:min(floor(numSamples/2),length(lowervals)));
            if length(testrangeL) < floor(numSamples/2)
                %add extra steps to have exactly numSamples/2 numbers in
                %the range
                l = length(testrangeL);
                testrangeL = [testrangeL,linspace(lbparam/5,lbparam,ceil(numSamples/2)-l)];
            end
        end
        testrange = [testrangeL,testrangeU];

        if lowersteps > 0 %only refine the closest to current limits
            disp('refine + lowersteps = re-estimate the n=lowersteps and n=uppersteps steps closest to current limits')
            testrange = [testrangeL(1:lowersteps),testrangeU(1:uppersteps)];
            testrange = repmat(testrange,1,ceil(numSamples/length(testrange)));% --> range is at least numSamples long
            testrange = testrange(1:numSamples);
            addRandomnessSmall=1; %only +-5% if addrandomness
        elseif increasebound
            disp('increasebound+refine = re-estimate the 6 steps closest to current limits, plus 2 extra steps')
            testrange = [ubparam+((max(highervals)-ubparam)/2), lbparam-((lbparam - min(lowervals))/2), testrangeL(1:3),testrangeU(1:3)];
            testrange = repmat(testrange,1,ceil(numSamples/length(testrange)));% --> range is at least numSamples long (ie if is 32 --> 4 repeats)
            testrange = testrange(1:numSamples);
            addRandomnessSmall=1; %only +-5% if addrandomness
        end

        %% Smooth: re-optimize steps where there is a sudden peak in the profile likelihood
    elseif smooth
        disp('SMOOOTH')
        %only smooth pvals with costs <300
        pvals = pvals(costs < 300);
        costs = costs(costs < 300);
        [peakCostvals,peakPvals] = findpeaks(costs,pvals,'MinPeakProminence',0.3);
        [peakCostvals2,peakPvals2] = findpeaks(-costs,pvals,'MinPeakProminence',0.3);
        [peakPvals,uind] = unique([peakPvals;peakPvals2]);
        peakCostvals =[peakCostvals;-peakCostvals2];
        peakCostvals = peakCostvals(uind);

        if isempty(peakPvals)
            disp('no peaks found!')
            break
        else
            testrange = peakPvals';
            testrange = repmat(testrange,1,ceil(numSamples/length(testrange)));% --> range is at least numSamples long (ie if is 32 --> 4 repeats)
            testrange = testrange(1:numSamples);
            %debug figure:
            % figure('Name','peaks')
            % plot(pvals,costs,'*')
            % hold on
            % plot(peakPvals,peakCostvals,'ro')
        end

        %% Increasebound: Load prev range to increase the range if needed
    elseif increasebound
        disp('Increase bounds')
        limit = bestcost+chi2inv(0.95,1)+chi2inv(0.95,1);
        withinrange = pvals(costs<=limit);
        lbparam = min(withinrange);
        ubparam = max(withinrange);

        newminval = NaN;
        newmaxval = NaN;
        if (lbparam-minval)/minval < 0.4 %minval close to edge of tested range
            newminval = max(1e-9, minval - minval*0.4);
        end
        if (maxval-ubparam)/maxval < 0.4  %maxval close to edge of tested range
            newmaxval = maxval + maxval*0.4;
        end

        step = 2*mean(diff(testrange));
        if lowersteps > 0
            newminval = min(pvals) - step*(lowersteps+1);
            if strcmp(paramname,'m2_LV')
                newminval = 1;
            end
        end

        if uppersteps > 0
            newmaxval = max(pvals) + step*(uppersteps+1);
        end

        if isnan(newminval) && ~isnan(newmaxval)
            testrange = linspace(maxval,newmaxval,numSamples);
        elseif isnan(newmaxval) && ~isnan(newminval)
            testrange = linspace(newminval,minval,numSamples);
        elseif ~isnan(newmaxval) && ~isnan(newminval)
            testrange1 = linspace(newminval,minval-step,ceil(numSamples/2));
            testrange2 = linspace(maxval+step,newmaxval,floor(numSamples/2));
            testrange = [testrange1,testrange2];
        else
            disp('no increase of bounds needed.')
            break
        end
    end

    %% Check that the testrange only includes reasonable values
    if sum(testrange < 100)>0
        disp('PL: removing parameter values <=0 from testrange, unphysiological value!')
        testrange(testrange <= 0) = [];
        testrange = repmat(testrange,1,ceil(numSamples/length(testrange)));% --> range is at least numSamples long
        testrange = testrange(1:numSamples); %make sure that the testrange still includes numSamples
    end

    if strcmp(paramname,'Cpvc') && sum(testrange > 100)>0
        disp('PL: Removing Cpvc > 100 from testrange, unphysiological value!')
        testrange(testrange > 100) = [];
        testrange = repmat(testrange,1,ceil(numSamples/length(testrange)));% --> range is at least numSamples long
        testrange = testrange(1:numSamples); %make sure that the testrange still includes numSamples
    end

    %% sort so that the parameter values closest to the optimal value are optimized first
    difffrombest = abs(testrange-optval);
    [s,sortind] = sort(difffrombest);
    testrange = testrange(sortind);

    %% set bounds
    datain = data;
    datain.meta.paramsToOptimize(indFixedParam) = [];% do not optimize the fixed param

    lb = log10(lbOrig);
    ub = log10(ubOrig);
    lbOrig(indFixedParam) = [];
    ubOrig(indFixedParam) = [];
    lb(indFixedParam) = [];
    ub(indFixedParam) = [];
    problem.x_L       = lb; % essOPT uses a problem structure where crucial information is specified
    problem.x_U       = ub;
    nParams = length(lb);

    %% find a bound for the parameter
    pvalue = testrange(loopNum);
    fprintf('-----Fixed %s to %0.4f-----------\n',paramname,pvalue)
    allparams(indFixedParam) = pvalue;
    %check if a paramset already optimized, and use as start guess
    savefolder = sprintf('%s/%s_%d',PLfolder,paramname,pvalue);
    bestcostAll = 1e29;
    for rep = 1:nRepeats
        loadfolders = dir([savefolder '/opt-*']);
        if isempty(loadfolders) %use best paramset as startguess
            mkdir(savefolder)
            optParam = bestParams;
            optParam(indFixedParam) = [];

            %increase bounds if optparam is outside the bounds
            if sum(optParam<lbOrig)>0 || sum(optParam>ubOrig)>0
                lbOrig = min(optParam,lbOrig);
                ubOrig = max(optParam,ubOrig);
                problem.x_L = log10(lbOrig); % essOPT uses a problem structure where crucial information is specified
                problem.x_U = log10(ubOrig);
            end
        else % load startguess from previous results
            [optParam,bestcostAll] = findBestParams(savefolder,0,nParams);
            if addRandomnessSmall && (addRandomness && rep < nRepeats/2)
                optParam = optParam.*( rand(1, length(optParam)).*0.1 + (1-(0.1/2)) );% +-5%
                optParam = max(optParam,lbOrig);
                optParam = min(optParam,ubOrig);
            elseif addRandomness && rep < nRepeats-10
                optParam = optParam.*(rand(1, length(optParam))./4+0.8750); % +-12.5%
                optParam = max(optParam,lbOrig);
                optParam = min(optParam,ubOrig);
            end
            disp('Loaded startguess')
        end
        if size(optParam) ~= size(lbOrig)
            optParam = optParam';
        end

        optParam =  log10(optParam);
        problem.x_0=optParam;
        warning('off','all') % AMICI prints error-messages for all simulations (if they fail) so this will fix that annoying orange text form appearing

        %% Solve with ESS
        optim_algorithm = 'ess';
        Results = MEIGO(problem,opts,optim_algorithm,constants,allparams,simOptions,datain,ind,dispErrors,doPlot); % Run the optimization

        %% Save results
        fitting_speed     = Results.time(end);
        bestcost           = Results.fbest;
        optParam          = Results.xbest;
        w = warning ('on','all'); % yay error messages again

        % Save the results.
        optParam = 10.^(optParam); % scale bestparam back to it's original size
        bounds.lb = 10.^lb;
        bounds.ub = 10.^ub;
        dateEnd = datestr(now,'yymmdd-HHMMSS');
        costChi2 = bestcost;
        if ~addRandomness && round(bestcost,2) == round(bestcostAll,2) && nRepeats>3 %if its not getting better, quit
            save(sprintf('%s/opt-PLbest.mat',savefolder) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
            disp('Cost is not improving, exiting PL for this parameter.')
            break;
        end
        if bestcost < bestcostAll
            bestcostAll = bestcost;
            save(sprintf('%s/opt-PLbest.mat',savefolder) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
        end
    end

end

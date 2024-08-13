function [] =  EstimatePL(loopNumber,dateStart,experimentNames,experimentRange,numberOfSamples,startN,paramToEstimate,doLower)
% Function to run profile likelihood for the bootstrapped data,
% starting at n=numSamples fixed steps, where each loopNum=1:numSamples optimizes one step.

%% Setup
s=rng('shuffle'); %do not note remove this line
d = datestr(now,'YYmmDD-HHMMSSFFF'); %creating different seed for each start
s=rng(s.Seed+loopNumber*sum(d)); %do not note remove this line

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

if nargin < 7
    paramNamesToEstimate = {'m2_LV'};
else
    paramNamesToEstimate = {paramToEstimate};
end


%% Load data
if ischar(experimentNames) || isstring(experimentNames)
    experimentNames = {experimentNames};
end

if strcmp(experimentNames{1}(1:3),'sim')% simulated data
    % Sample data
    n = mod(loopNumber,numberOfSamples)+1+startN;
    expName = experimentNames{1};
    [data,~] = loadSampledMeasurementError(n,expName);
    data.(expName) = data; %create struct
    simulatedData = 1;
    disp(['----------------------- Exp: ' num2str(n) ' ---------------------------------'])
else % "measured" or other data
    simulatedData = 0;
    % Choose which experiment to simulate based on the range of experiments and which loop number you are in
    % Load the data to get all experiment names
    load('data.mat','data')
    % select experiments
    [~,experimentNames] = extractData(experimentNames,experimentRange,data);

    num = mod(loopNumber,length(experimentNames))+1;
    expName = experimentNames{num};
    disp(['----------------------- Exp: ' expName ' ---------------------------------'])
end


if nargin < 8
    % do lower or upper based on the loopnumber. 
    % Assumes that loopnumber is 1:numberOfSamples*2
    if loopNumber > numberOfSamples
        doLower = 0;
    else
        doLower = 1;
    end
end

%% Load parameters and data
resultsfolder = fullfile(basefolder,'Parameters'); 
[~,data,~,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues,~,~,~,~,~] = setup_simulations({expName},data,resultsfolder,1);

% settings for the cost function:
doPlot = 0;
dispErrors = 0;


%% number of optimizations for each parameter value
if strcmp(data.meta.simulationFunction,'simulate_avatar_correction')
    nRepeats = 5; %is much slower to simulate, thus fewer optimizations per step
else
    nRepeats = 7;
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
if simulatedData
    folders = dir(fullfile(resultsfolder,'ESS',['E_' expName num2str(n) '_*']));
else
    folders = dir(fullfile(resultsfolder,'ESS',['E_' expName '*']));
end

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

%%
thresholdChi2 = chi2inv(0.95,1); %if only for this parameter

%% find folder to save results in for this data
if simulatedData
    PLfolder = fullfile(basefolder,'Parameters','PL',['E_' expName num2str(n)]);
else
    PLfolder = fullfile(basefolder,'Parameters','PL', ['E_' expName]);
end

%% Run optimization
format long
format compact
warning('off','all') % AMICI prints error-messages for all simulations (if they fail) so this will fix that annoying orange text form appearing

for n = 1:length(paramNamesToEstimate)
    paramname = paramNamesToEstimate{n};
    indFixedParam = strcmp(paramname,paramNames);
    fprintf('Starting paramestimation of %s for %s \n',paramname,expName)

    %% Set test ranges
    minval = lbOrig(indFixedParam)/5;
    optval = bestParams(indFixedParam);
    maxval = ubOrig(indFixedParam)*5;
    step  = (maxval-minval)/50;
    if doLower
        testrange = optval:-step:minval;
        loweruppername = 'Lower';
    else
        testrange = optval+step:step:maxval;
        loweruppername = 'Upper';
    end

    %% Special case for onset_LV
    if strcmp(paramname,'onset_LV')
        disp('Special case: fixed testrange for onset_LV')
        minval = 1.2;
        maxval = 1.6;
        step  = (maxval-minval)/10;
        if doLower
            testrange = optval:-step:minval;
        else
            testrange = optval+step:step:maxval;
        end
    end

    %% check if PL was previously done, and use those steps in that case
    prevfolders = dir(sprintf('%s/%s_*',PLfolder,paramname));
    prevtestrange = nan(1,length(prevfolders));
    for i = 1:length(prevfolders)
        savefolder = fullfile(prevfolders(i).folder,prevfolders(i).name);
        if ~isempty(dir([savefolder '/opt-*']))
            val = split(prevfolders(i).name,'_');
            val = str2double(val{end});
            prevtestrange(i) = val;
        end
    end
    prevtestrange = sort(prevtestrange);
    prevtestrange(isnan(prevtestrange)) = [];
    if ~isempty(prevtestrange)
        %replace the steps in the testrange with the prev values to be
        %able to load a good start guess

        %find the previously  optimized steps that are within the current
        %wanted testrange
        withinrange = (prevtestrange >= min(testrange)) & (prevtestrange <= max(testrange));
        rangetoadd = prevtestrange(withinrange);

        %remove steps that are within the same range as the already optimized before
        testrange(testrange <= max(rangetoadd) & testrange >= min(rangetoadd)) = [];

        %add the previously optimized steps instead
        if doLower
            testrange = [testrange,rangetoadd,optval]; %make sure that the optval is included and it goes from larger to smaller values away from optimum
            testrange=unique(testrange);
            testrange = sort(testrange,'descend');
        else
            testrange = sort([testrange,rangetoadd]);
        end
    end

    %% set bounds
    datain = data;
    datain.meta.paramsToOptimize(indFixedParam) = [];% do not optimize the fixedx param

    lb = log10(lbOrig);
    ub = log10(ubOrig);
    lb(indFixedParam) = [];
    ub(indFixedParam) = [];
    problem.x_L       = lb; % essOPT uses a problem structure where crucial information is specified
    problem.x_U       = ub;
    nParams = length(lb);
    foundbound = 0;

    %% find a bound for the parameter
    allcosts = nan(length(testrange)+1,2);

    i = 1;
    for pvalue = testrange
        fprintf('-----Fixed %s to %0.4f-----------\n',paramname,pvalue)
        allparams(indFixedParam) = pvalue;
        %check if a paramset already optimized, and use as start guess
        savefolder = sprintf('%s/%s_%d',PLfolder,paramname,pvalue);
        bestcostAll = 1e29;
        for rep = 1:nRepeats
            loadfolders = dir([savefolder '/opt-*']);
            if isempty(loadfolders) %use best paramset as startguess
                mkdir(savefolder)
                if i == 1 %start at optimum
                    optParam = bestParams;
                    optParam(indFixedParam) = [];
                    %othewise start at optimum for prev. pvalue (optParam)
                end
            else % load startguess from previous results 
                [optParam,optcost] = findBestParams(savefolder,0,nParams);
                %optParam = optParam.*(rand(1, length(optParam))./8+0.9375);%./4+0.8750
                disp('Loaded startguess')
            end
            if size(optParam) ~= size(lb)
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
            save(sprintf('%s/opt-PL(%.3f)-%i.mat',savefolder,bestcost,loopNumber) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
            if round(bestcost,3) == round(bestcostAll,3) %if its not getting better, quit
                bestcostAll = bestcost;
                break;
            end
            if bestcost < bestcostAll
                bestcostAll = bestcost;
            end
        end
        allcosts(i,2) = bestcostAll;
        allcosts(i,1) = pvalue;
        save(sprintf('%s/allcosts%s_%s_%i.mat',PLfolder,loweruppername,paramname,loopNumber),'allcosts') %save each round in case it craches  

        %update the minimum and the limit if a better cost was found
        if bestcostAll < optcostChi2
            optcostChi2 = bestcostAll;
        end

        %quit if found a cost over the threshold for this parameter (ie
        %a upper or lower bound)  (+ one degrees of freedom to make sure)
        if foundbound && bestcostAll > optcostChi2 + thresholdChi2+chi2inv(0.95,1)
            disp('Paramer bound found, stopping.')
            i= i+1;
            break;
        elseif bestcostAll > optcostChi2 + thresholdChi2+chi2inv(0.95,1)
            foundbound = 1;
        end
        i= i+1;
    end
    save(sprintf('%s/allcosts%s_%s_%i.mat',PLfolder,loweruppername,paramname,loopNumber),'allcosts')

    %% do one more optimization at the in-between step if found the bound
    if foundbound && i>3 
        pvalue = (testrange(i-2) + testrange(i-3))/2;%take half a step back %i ends at one more than the last
        fprintf('-----Fixed %s to %0.4f-----------\n',paramname,pvalue)
        allparams(indFixedParam) = pvalue;
        %check if a paramset already optimized, and use as start guess
        savefolder = sprintf('%s/%s_%d',PLfolder,paramname,pvalue);
        bestcostAll = 1e29;
        loadfolders = dir([savefolder '/opt-*']);
        if isempty(loadfolders) %use best paramset as startguess
                mkdir(savefolder)
        end
        for rep = 1:nRepeats
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
            optParam = 10.^(optParam); % scale bestparam back to it's original size
            bounds.lb = 10.^lb;
            bounds.ub = 10.^ub;
            dateEnd = datestr(now,'yymmdd-HHMMSS');
            costChi2 = bestcost;
            save(sprintf('%s/opt-PL(%.3f)-%i.mat',savefolder,bestcost,loopNumber) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
            if round(bestcost,3) == round(bestcostAll,3) %if its not getting better, quit
                bestcostAll = bestcost;
                break;
            end
            if bestcost < bestcostAll
                bestcostAll = bestcost;
            end
        end
        allcosts = [allcosts;[pvalue,bestcostAll]];
        save(sprintf('%s/allcosts%s_%s_%i.mat',PLfolder,loweruppername,paramname,loopNumber),'allcosts')
    end
    save(sprintf('%s/allcosts%s_%s_%i_finished.mat',PLfolder,loweruppername,paramname,loopNumber),'allcosts')

end

end

function [] =  EstimatePL_realdata_classic(loopNum,dateStart,dataName,continueFromLimit,paramToEstimate,addRandomness)
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


if strcmp(paramToEstimate,'') || isempty(paramToEstimate)
    if strcmp(dataName,'dataP33')
        paramNamesToEstimate ={'Rao','m2_LV','m2_LV','Emax_LA','Emax_LA','k_syst_LV','k_syst_LV','Caa','Caa','Rao'};
    else
        paramNamesToEstimate ={'Rao','m2_LV','m2_LV','Cpvc','Cpvc','Emax_LA','Emax_LA','k_syst_LV','k_syst_LV','Caa','Caa','Rao'};
    end
else
    paramNamesToEstimate = {paramToEstimate};
end


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
if addRandomness
    nRepeats = 20;   %extra in the end to improve the best found (last half of reps are not using random startguess)
else
    nRepeats = 10;
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



%% Set bounds (separate for each subject, since paramdata is individual)
[lbOrig,ubOrig] = loadParamBounds(ind,data.parameters);
lbOrig = lbOrig(data.meta.paramsToOptimize);
ubOrig = ubOrig(data.meta.paramsToOptimize);

allparams = origParamvalues;

%% find folder to save results in for this data
PLfolder = fullfile(basefolder,'Parameters','PL',['E_' dataName]);


format long
format compact
warning('off','all') % AMICI prints error-messages for all simulations (if they fail) so this will fix that annoying orange text form appearing

%% find parameter to estimate
if length(paramNamesToEstimate)>1
    n = mod(loopNum, length(paramNamesToEstimate))+1;
else
    n=1;
end

paramname = paramNamesToEstimate{n};
indFixedParam = strcmp(paramname,paramNames);
fprintf('Starting paramestimation of %s for %s \n',paramname,dataName)

%% set bounds
datain = data;
datain.meta.paramsToOptimize(indFixedParam) = [];% do not optimize the fixedx param

lb = log10(lbOrig);
ub = log10(ubOrig);
lbOrig(indFixedParam) = [];
ubOrig(indFixedParam) = [];
lb(indFixedParam) = [];
ub(indFixedParam) = [];
problem.x_L       = lb; % essOPT uses a problem structure where crucial information is specified
problem.x_U       = ub;
nParams = length(lb);

%% find prev optimized range
% load all PLs for this experiment
loadParameters=1;% loading best parameters for the subject from all saved parameter values
[foundPL,allPLs,~,~,bestParams,bestSubjectCost] = findPLparams(dataName,resultsfolder,length(paramNames),loadParameters,{paramname},paramNames);
if ~foundPL
    disp('EstimatePL_realdata_classic: did not find PL params, EXITING')
    pvals = [];
elseif length(allPLs.(paramname)(:,1))<2
    disp('EstimatePL_realdata_classic: no PL params, using default trange')
        optval = bestParams(indFixedParam);
    if length(bestParams) > 28
        bestParams = bestParams(1:28);
    end
    bestParams(indFixedParam) = [];
    step = optval*0.1;
    pvals = optval-(step*40):step:optval+(step*40);
    pvals(pvals<=0) = [];

    nRepeats = 40;
    disp('EstimatePL_realdata_classic: increased nRepeats to due to the new PL')
else
    optval = bestParams(indFixedParam);
    if length(bestParams) > 28
        bestParams = bestParams(1:28);
    end
    bestParams(indFixedParam) = [];

    %sort
    [~,sind]=sort(allPLs.(paramname)(:,1));
    costs = allPLs.(paramname)(sind,2);
    pvals = allPLs.(paramname)(sind,1);
    % [bestcost,costind] = min(costs);
    % optval = pvals(costind);

    if continueFromLimit %do not start at optimum, start further away close to the limits of the confidence interval
        disp('EstimatePL_realdata_classic: continueFromLimit')
        limit = bestSubjectCost+chi2inv(0.95,1);
        maxval = max(pvals(costs<limit));
        minval = min(pvals(costs<limit));
        pvals  = pvals(costs >= limit);
        pvals = [pvals;maxval;minval];
        pvals = sort(pvals);
    end

    if ~strcmp('Cpvc',paramname) %length(pvals)<40
        disp('EstimatePL_realdata_classic: adding +-10 extra steps with 4* the mean step length')
        step = mean(diff(pvals));
        try
            pvals = [(pvals(1)-(step*40)):step*4:pvals(1), pvals, pvals(end):step*4:(pvals(end)+step*40)];
        catch
            pvals = [(pvals(1)-(step*40)):step*4:pvals(1), pvals', pvals(end):step*4:(pvals(end)+step*40)]';
        end
        pvals = unique(pvals);
        pvals(pvals<=0) = [];
    end

end

if rem(loopNum, 2) == 0 %even
    testrange = pvals(pvals <= optval);%lower
    if length(testrange)<3
        testrange = [testrange;min(testrange)*0.8;min(testrange)*0.5;min(testrange)*0.1;min(testrange)*0.01];
    end
    testrange = sort(testrange, 'descend');
    dolower = 1;
    disp('EstimatePL_realdata_classic: doLower')
else %odd
    testrange = pvals(pvals >= optval);%higher
    if length(testrange)<2 && ~strcmp('Cpvc',paramname)
        testrange = [testrange;max(testrange)*1.2;max(testrange)*1.5;max(testrange)*2;max(testrange)*5;max(testrange)*10];
    end
    testrange = sort(testrange, 'ascend');
    dolower = 0;
    disp('EstimatePL_realdata_classic: doUpper')
end

% remove unphysiological values
if strcmp('Cpvc',paramname)  && sum(testrange>100)>0
    testrange = testrange(testrange<=100);
    testrange = [testrange,100];
    disp('EstimatePL_realdata_classic: removing Cpvc > 100')
elseif strcmp('k_syst_LV',paramname) && sum(testrange>200)>0
    testrange = testrange(testrange<=200);
    testrange = [testrange,200];
    disp('EstimatePL_realdata_classic: removing k_syst_LV > 200')
end

%% Run optimization
for t = 1:length(testrange)
    %optimize the step using the closest step as the first startguess
    pvalue = testrange(t);
    fprintf('-----Fixed %s to %0.4f-----------\n',paramname,pvalue)
    allparams(indFixedParam) = pvalue;
    %check if a paramset already optimized, and use as start guess
    savefolder = sprintf('%s/%s_%d',PLfolder,paramname,pvalue);
    bestcostAll = 1e29;
    for rep = 1:nRepeats
        loadfolders = dir([savefolder '/opt-*']);
        if rep == 1 && continueFromLimit && ~isempty(loadfolders)% load startguess from previous results
            [optParam,bestcostAll] = findBestParams(savefolder,0,nParams);
            if size(optParam) ~= size(lbOrig)
                optParam = optParam';
            end
            disp('Loaded startguess')
        elseif rep == 1 || isempty(loadfolders) %always use the best paramset from the previous step as the first startguess
            disp('Using best param as startguess')
            mkdir(savefolder)
            optParam = bestParams;
            if size(optParam) ~= size(lbOrig)
                optParam = optParam';
            end
            %increase bounds if optparam is outside the bounds
            if sum(optParam<lbOrig)>0 || sum(optParam>ubOrig)>0
                lbOrig = min(optParam,lbOrig);
                ubOrig = max(optParam,ubOrig);
                problem.x_L = log10(lbOrig);
                problem.x_U = log10(ubOrig);
            end
        elseif rep <4 %continue with the previous reps opt for rep 2 and 3
            if addRandomness && rep < nRepeats/2
                optParam = optParam.*( rand(1, length(optParam)).*0.1 + (1-(0.1/2)) );% +-5%
                if size(optParam) ~= size(lbOrig)
                    optParam = optParam';
                end
                optParam = max(optParam,lbOrig);
                optParam = min(optParam,ubOrig);
            end
            disp('Continuing on previous repetition as startguess')
        else % load startguess from previous results
            [optParam,bestcostAll] = findBestParams(savefolder,0,nParams);
            if size(optParam) ~= size(lbOrig)
                optParam = optParam';
            end
            if addRandomness && rep < nRepeats/2
                optParam = optParam.*( rand(1, length(optParam)).*0.1 + (1-(0.1/2)) );% +-5%
                % optParam = optParam.*(rand(1, length(optParam))./4+0.8750); % +-12.5%
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

        % Solve with ESS
        optim_algorithm = 'ess';
        Results = MEIGO(problem,opts,optim_algorithm,constants,allparams,simOptions,datain,ind,dispErrors,doPlot); % Run the optimization

        % Save results
        fitting_speed     = Results.time(end);
        bestcost           = Results.fbest;
        optParam          = Results.xbest;
        w = warning ('on','all'); % yay error messages again

        optParam = 10.^(optParam); % scale bestparam back to it's original size
        bounds.lb = 10.^lb;
        bounds.ub = 10.^ub;
        dateEnd = datestr(now,'yymmdd-HHMMSS');
        costChi2 = bestcost;

        %check that no better costs have been found elsewhere
        load(sprintf('%s/opt-PLbest.mat',savefolder) ,'bestcost')
        if bestcost < bestcostAll
            bestcostAll=bestcost;
            disp('loaded bestcost from another optimization arm')
        end
        bestcost=costChi2;%back to bestcost from this optimization

        if bestcost < bestcostAll
            bestcostAll = bestcost;
            save(sprintf('%s/opt-PLbest.mat',savefolder) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
            fprintf('saved best params in %s/opt-PLbest.mat\n',savefolder)
        end

        save(sprintf('%s/opt-PL(%.3f)-%i.mat',savefolder,bestcost,loopNum) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
        if ~addRandomness && round(bestcost,2) == round(bestcostAll,2) && nRepeats>2 %if its not getting better, quit
            %save(sprintf('%s/opt-PLbest.mat',savefolder) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
            disp('Cost is not improving, exiting PL for this parameter.')
            break;
        elseif addRandomness  && round(bestcost,2) == round(bestcostAll,2) && nRepeats>(nRepeats/2) %if its not getting better, quit
            %save(sprintf('%s/opt-PLbest.mat',savefolder) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
            disp('Cost is not improving, exiting PL for this parameter.')
            break;
        end

    end

    %use previous best parameters as startguess for the next step
    bestParams = optParam;

    if bestcostAll < bestSubjectCost
        bestSubjectCost = bestcostAll;
    end

    limit = bestSubjectCost+chi2inv(0.95,1);
    if bestcostAll > limit
        disp('EstimatePL_realdata_classic: bestcostAll > limit, exiting from testrange.')
        % continue to the next step if you are still below the limit,
        % otherwise stop and add another step in between and optimize that one
        break
    end

end


%% Final round
if t > 1
    diffpval = abs(pvalue-testrange(t-1));
else
    diffpval = min(diff(pvals));
end
if bestcostAll > limit
    %add another step in between this step and the previous step, and optimize that one
    if dolower
        pvalue = pvalue+(diffpval/2);
    else
        pvalue = pvalue-(diffpval/2);
    end
else
    % add a larger/smaller step
    if dolower
        pvalue = pvalue-diffpval;
    else
        pvalue = pvalue+diffpval;
    end
end


fprintf('-----Final round: Fixed %s to %0.4f-----------\n',paramname,pvalue)
allparams(indFixedParam) = pvalue;
savefolder = sprintf('%s/%s_%d',PLfolder,paramname,pvalue);
bestcostAll = 1e29;
for rep = 1:nRepeats
    loadfolders = dir([savefolder '/opt-*']);
    if rep == 1 || isempty(loadfolders) %always use the best paramset from the previous step as the first startguess
        mkdir(savefolder)
        optParam = bestParams;
        %increase bounds if optparam is outside the bounds
        if sum(optParam<lbOrig)>0 || sum(optParam>ubOrig)>0
            lbOrig = min(optParam,lbOrig);
            ubOrig = max(optParam,ubOrig);
            problem.x_L = log10(lbOrig);
            problem.x_U = log10(ubOrig);
        end
    else % load startguess from previous results
        [optParam,bestcostAll] = findBestParams(savefolder,0,nParams);
        if addRandomness && rep < nRepeats/2
            optParam = optParam.*( rand(1, length(optParam)).*0.1 + (1-(0.1/2)) );% +-5%
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

    % Solve with ESS
    optim_algorithm = 'ess';
    Results = MEIGO(problem,opts,optim_algorithm,constants,allparams,simOptions,datain,ind,dispErrors,doPlot); % Run the optimization

    % Save results
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


    %check that no better costs have been found elsewhere
    load(sprintf('%s/opt-PLbest.mat',savefolder) ,'bestcost')
    if bestcost < bestcostAll
        bestcostAll=bestcost;
        disp('loaded bestcost from another optimization arm')
    end
    bestcost=costChi2;%back to bestcost from this optimization

    if bestcost < bestcostAll
        bestcostAll = bestcost;
        save(sprintf('%s/opt-PLbest.mat',savefolder) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
        fprintf('saved best params in %s/opt-PLbest.mat\n',savefolder)
    end

    % save(sprintf('%s/opt-PL(%.3f)-%i.mat',savefolder,bestcost,loopNum) ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNamesToEstimate','dateStart','dateEnd','fitting_speed','s','Results','pvalue')
    if ~addRandomness && round(bestcost,2) == round(bestcostAll,2) && nRepeats>3 %if its not getting better, quit
        disp('Cost is not improving, exiting PL for this parameter.')
        break;
    end

end




end

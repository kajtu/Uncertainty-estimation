function [] =  EstimatePL_refine(loopNumber,dateStart,experimentNames,experimentRange,numberOfSamples,startN,paramToEstimate)
% Function to run profile likelihood for the bootstrapped data,
% starting at n=numSamples fixed steps, where each loopNum=1:numSamples optimizes one step.
% The _refine function specifically adds more steps close to the true value

%% Setup
s=rng('shuffle'); %do not note remove this line
d = datestr(now,'YYmmDD-HHMMSSFFF'); %creating different seed for each start
s=rng(s.Seed+loopNumber*sum(d)); %do not note remove this line

basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');

%add all project folders to the matlab path
addpath(genpath(basefolder))

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
    [data,errorn] = loadSampledMeasurementError(n,expName);
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


%% Load parameters and data
resultsfolder = fullfile(basefolder,'Parameters'); 
[~,data,~,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues,~,~,~,~,~] = setup_simulations({expName},data,resultsfolder,1);

% settings for the cost function:
doPlot = 0;
dispErrors = 0;

%% Load true data
load('dataSimulated.mat','dataSimulated','simulatedDataTrue')
trueParams = simulatedDataTrue.allParameters;


%% ESS options
% MEIGO OPTIONS I (COMMON TO ALL SOLVERS):
opts.ndiverse     = 'auto'; 
if strcmp(data.meta.simulationFunction,'simulate_avatar_correction')
    opts.maxtime      = 2*600; % MAX-Time of optmization, i.e how long the optimization will last %Maximum CPU time in seconds (Default 60)
    nRepeats = 3;
else
    opts.maxtime      = 600; % MAX-Time of optmization, i.e how long the optimization will last %Maximum CPU time in seconds (Default 60)
    nRepeats = 5;
end
opts.maxeval      = 1e8; % max number of evals, i.e cost function calls (Default 1000)
opts.log_var      = [];  %skip this

opts.local.solver = 'dhc'; 
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

    %% check if refinement already done
    paramfilesPrevious = dir(sprintf('%s/allcosts*_%s_*ref_finished3*',PLfolder,paramname));

    if ~isempty(paramfilesPrevious)
        disp('this parameter is already refined! Exiting PL_refine...')
        break
    end

    %% load previous results to improve on
    %check if there are finished lower & upper estimates that needs to be
    %combined
    foundPL = false;

    paramfiles = dir(sprintf('%s/allcosts*_%s*',PLfolder,paramname));
    combinedcosts = [];
    for f = 1:length(paramfiles)
        load(fullfile(paramfiles(f).folder,paramfiles(f).name),'allcosts')
        combinedcosts = [combinedcosts;allcosts];
    end
    if ~isempty(combinedcosts)
        [~,indsu] = unique(combinedcosts(:,1));
        combinedcostsUnique=combinedcosts(indsu,:);
        combinedcostsUnique = combinedcostsUnique(~isnan(combinedcostsUnique(:,1)),:);
        for j = 1:length(combinedcostsUnique(:,1)) %make sure that the lowest cost for each pvalue is found
            pval = combinedcostsUnique(j,1);
            pcosts = combinedcosts(combinedcosts(:,1) == pval,2);
            combinedcostsUnique(j,2) = min(pcosts);
        end
        foundPL = true;
        combinedcosts = combinedcostsUnique;
    end


    if foundPL
        prevPvals = combinedcosts(:,1);
        prevCosts = combinedcosts(:,2);
    else
        disp('OBS no prev PL params found!!!')
        break
    end

    %% compare vals to trueval to refine
    trueval = trueParams(indFixedParam);

    if strcmp(paramname,'mvCorr')
        trueval = errorn.mv.systematic+40;
    end

    if min(prevPvals) > trueval
        %prev pl not low enough, optimize around trueval
        diff = min(prevPvals)-trueval;
        if trueval-(diff/2) <= 0
            testrange =  [trueval-(diff/6), trueval-(diff/3)];
            testrange(testrange <= 0) = trueval-0.0001;
        elseif trueval-diff <= 0
            testrange =  [trueval-(diff/4), trueval-(diff/2)];
        else
            testrange =  [trueval-(diff/2), trueval-diff];%no values below 0
        end
        closestPval = min(prevPvals);
        loweruppername = 'Lower';
    elseif max(prevPvals) < trueval
        %prev pl not high enough, optimize around trueval
        diff = trueval - max(prevPvals);
        testrange =  [trueval+(diff/2), trueval+diff];
        closestPval = max(prevPvals);
        loweruppername = 'Upper';
    else
        %optimize in between the vals
        aboveTrueval = prevPvals(prevPvals>trueval);
        belowTrueval = prevPvals(prevPvals<trueval);
        %prev vals closest to trueval (above)
        [mindist,minind1] = min(abs(aboveTrueval-trueval));
        minvalAbove = aboveTrueval(minind1);
        %prev vals closest to trueval (below)
        [mindist,minind2] = min(abs(belowTrueval - trueval));
        minvalBelow = belowTrueval(minind2);
        
        step = abs(minvalBelow-minvalAbove)/3;
        testrange =  [minvalBelow+step, minvalBelow+(step*2)];
        closestPval = minvalBelow;
        loweruppername = 'Lower';
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

    %% Load prev. results for the closest param value
    prevfolder = sprintf('%s/%s_%d',PLfolder,paramname,closestPval);
    [optParam,optcost] = findBestParams(prevfolder,0,nParams);

    %% find a bound for the parameter
    allcostsRef = nan(length(testrange),2);
    i = 1;
    for pvalue = testrange
        fprintf('-----Fixed %s to %0.4f-----------\n',paramname,pvalue)
        allparams(indFixedParam) = pvalue;
        %check if a paramset already optimized, and use as start guess
        savefolder = sprintf('%s/%s_%d',PLfolder,paramname,pvalue);
        bestcostAll = 1e29;
        stopAfterOneOpt = 0;
        for rep = 1:nRepeats
            loadfolders = dir([savefolder '/opt-*']);
            if isempty(loadfolders) %use best prev sparamset as startguess
                mkdir(savefolder)
            else % load startguess from prev results in this run
                [optParam,optcost,~,~,~,~,~,bestcostAll] = findBestParams(savefolder,0,nParams);%[bestparams,bestcostAll,allparams,bestconstants,meanparams,medianparams,allCostsChi2,bestcostAllchi2]
                disp('Loaded startguess')

                if rep == 1
                    disp('paramvalue already opt, only doing one more opt for this pvalue')
                    stopAfterOneOpt = 1;
                end
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
            if round(bestcost,3) == round(bestcostAll,3) || stopAfterOneOpt %if its not getting better
                bestcostAll = bestcost;
                break;
            end
            if bestcost < bestcostAll
                bestcostAll = bestcost;
            end
        end
        allcostsRef(i,2) = bestcostAll;
        allcostsRef(i,1) = pvalue;
        allcosts = [allcostsRef;combinedcosts];
        i= i+1;
    end
    allcosts = [allcostsRef;combinedcosts];

    save(sprintf('%s/allcosts%s_%s_%i_ref_finished3.mat',PLfolder,loweruppername,paramname,loopNumber),'allcosts')

end

end

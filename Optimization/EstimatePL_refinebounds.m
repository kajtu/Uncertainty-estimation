function [] =  EstimatePL_refinebounds(loopNumber,dateStart,experimentNames,experimentRange,numberOfSamples,startN,paramToEstimate)
% Function to run profile likelihood for the bootstrapped data,
% starting at n=numSamples fixed steps, where each loopNum=1:numSamples optimizes one step.
% The _refinebounds function specifically adds more steps close to the chi2
% limits.

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
nRepeats = 7;

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

% do lower or upper based on the loopnumber. 
% Assumes that loopnumber is 1:numberOfSamples*2
if loopNumber > numberOfSamples
    doLower = 0;
else
    doLower = 1;
end

%% Load parameters and data
resultsfolder = fullfile(basefolder,'Parameters'); 
[~,data,~,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues,~,~,~,~,~] = setup_simulations({expName},data,resultsfolder,1);

% settings for the cost function:
doPlot = 0;
dispErrors = 0;


%% ESS options
% MEIGO OPTIONS I (COMMON TO ALL SOLVERS):
opts.ndiverse     = 'auto'; 
if strcmp(data.meta.simulationFunction,'simulate_avatar_correction')
    opts.maxtime      = 2*600; % MAX-Time of optmization, i.e how long the optimization will last %Maximum CPU time in seconds (Default 60)
else
    opts.maxtime      = 600; % MAX-Time of optmization, i.e how long the optimization will last %Maximum CPU time in seconds (Default 60)
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
    [bestParamsESS,optcost,~,loadedconstants,~,~,~,optcostChi2] = findBestParams(loadfolderName,0,length(paramNames));
    fprintf('Loaded paramvalues with cost %0.2f\n',optcostChi2)
end

%% Set bounds (separate for each subject, since paramdata is individual)
[lbOrig,ubOrig] = loadParamBounds(ind,data.parameters);
lbOrig = lbOrig(data.meta.paramsToOptimize);
ubOrig = ubOrig(data.meta.paramsToOptimize);

allparams = origParamvalues;

%% find folder to save results in for this data
if simulatedData
    folder = fullfile(basefolder,'Parameters','PL',['E_' expName num2str(n)]);
else
    folder = fullfile(basefolder,'Parameters','PL', ['E_' expName]);
end

%% Run optimization
format long
format compact
warning('off','all') % AMICI prints error-messages for all simulations (if they fail) so this will fix that annoying orange text form appearing

for n = 1:length(paramNamesToEstimate)
    paramname = paramNamesToEstimate{n};
    indFixedParam = strcmp(paramname,paramNames);
    fprintf('Starting paramestimation of %s for %s \n',paramname,expName)

    %% FIND PREVIOUS PL SAMPLING FOR THIS PARAMETER
    [foundPL,allPL] = findPLparams([expName num2str(n)],fullfile(basefolder,'Parameters'),length(paramNames),0,paramNames);
    if ~foundPL
        break
    else
        % sort
        [~,sind]=sort(allPL.(paramname)(:,1));
        costs = allPL.(paramname)(sind,2);
        pvals = allPL.(paramname)(sind,1);
        allcostsLoaded = allPL.(paramname)(sind,:);
        % calculate limit
        [mincost, midind] = min(costs);
        limit = mincost+chi2inv(0.95,1);
        % find bounds
        ub = max(pvals(costs <= limit));
        lb = min(pvals(costs <= limit));
        % define new test ranges
        if doLower && lb == min(pvals)  || ~doLower && ub == max(pvals)
            disp('The parameter is unidentifiable for the testrange')
            break
        elseif doLower
            pvalUnder = lb;
            pvalOver = pvals(find(pvals == pvalUnder)-1);
            minVal = pvalOver;
            maxVal = pvalUnder;
        else 
            pvalUnder = ub;
            pvalOver = pvals(find(pvals == pvalUnder)+1);
            minVal = pvalUnder;
            maxVal = pvalOver;
        end

        % calculate test range
        dc =  costs(pvals == maxVal)  - costs(pvals == minVal);
        dp = maxVal - minVal;
        crange = [limit-1:0.1:limit+1];
        testrange = minVal + ((crange - costs(pvals == minVal))/(dc/dp));
    end

    if doLower
        loweruppername = 'Lower';
        lbfolder = sprintf('%s/%s_%d',folder,paramname,lb);
        try
            [bestParams,optcost] = findBestParams(lbfolder,0,nParams);
        catch
            disp('could not load pl params for lb')
            bestParams = bestParamsESS;
        end
    else
        loweruppername = 'Upper';
        ubfolder = sprintf('%s/%s_%d',folder,paramname,ub);
        try
            [bestParams,optcost] = findBestParams(ubfolder,0,nParams);
        catch
            disp('could not load pl params for ub')
            bestParams = bestParamsESS;
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

    %% find a bound for the parameter
    save(sprintf('%s%sEstimatePL_refine_start%s_%s-%i-%s.mat',folder,filesep,loweruppername,paramname,loopNumber,dateStart) ,'dateStart')
    allcostsRef = nan(length(testrange),2);
    i = 1;
    for pvalue = testrange
        fprintf('-----Fixed %s to %0.4f-----------\n',paramname,pvalue)
        allparams(indFixedParam) = pvalue;
        %check if a paramset already optimized, and use as start guess
        savefolder = sprintf('%s/%s_%d',folder,paramname,pvalue);
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
            else % load startguess from prev results in this run
                [optParam,optcost] = findBestParams(savefolder,0,nParams);
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
            if round(bestcost,3) == round(bestcostAll,3) %if its not getting better
                bestcostAll = bestcost;
                break;
            end
            if bestcost < bestcostAll
                bestcostAll = bestcost;
            end
        end
        allcostsRef(i,2) = bestcostAll;
        allcostsRef(i,1) = pvalue;
        allcosts = [allcostsRef;allcostsLoaded];
        save(sprintf('%s/allcosts%s_%s_%i.mat',folder,loweruppername,paramname,loopNumber),'allcosts') %save each round in case it craches 
        i= i+1;
    end
    allcosts = [allcostsRef;allcostsLoaded];
    save(sprintf('%s/allcosts%s_%s_%i_finished.mat',folder,loweruppername,paramname,loopNumber),'allcosts')

end

end

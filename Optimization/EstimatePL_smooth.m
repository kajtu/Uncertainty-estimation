function [] =  EstimatePL_smooth(loopNumber,dateStart,experimentNames,experimentRange,numberOfSamples,startN,paramToEstimate)
% Function to run profile likelihood for the bootstrapped data,
% starting at n=numSamples fixed steps, where each loopNum=1:numSamples optimizes one step.
% The _smooth function specifically re-optimizes all steps that are peaks
% in the profile likelihood for each biomarker

%% Setup
s=rng('shuffle'); %do not note remove this line
d = datestr(now,'YYmmDD-HHMMSSFFF'); %creating different seed for each start
s=rng(s.Seed+loopNumber*sum(d)); %do not note remove this line

disp(['inside EstimatePL_smooth! ' num2str(loopNumber)])

basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');

% Create log file
% fid=fopen(fullfile(basefolder,'Log',['log_estPLsmooth_' dateStart '_' paramToEstimate num2str(loopNumber) '.txt']),'a');
% fprintf(fid, ['inside EstimatePL_smooth!' num2str(loopNumber) '\n']);

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
    expNameunique = [expName num2str(n)];

else % "measured" or other data
    simulatedData = 0;
    % Choose which experiment to simulate based on the range of experiments and which loop number you are in
    % Load the data to get all experiment names
    load('data.mat','data')
    % select experiments
    [~,experimentNames] = extractData(experimentNames,experimentRange,data);

    num = mod(loopNumber,length(experimentNames))+1;
    expName = experimentNames{num};
    expNameunique= expName;
    disp(['----------------------- Exp: ' expName ' ---------------------------------'])
end



%% Load parameters and data
resultsfolder = fullfile(basefolder,'Parameters'); 
[~,data,~,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues,~,~,~,~,~] = setup_simulations({expName},data,resultsfolder,1);

% settings for the cost function:
doPlot = 0;
dispErrors = 0;

%% number of optimizations for each parameter value
if strcmp(data.meta.simulationFunction,'simulate_avatar_correction')
    nRepeats = 2; %is much slower to simulate, thus fewer optimizations per step
else
    nRepeats = 4;
end

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

    loadParameters = 1;%load parameters from folders instead of combinedcosts summary files
    [foundPL,allPLs,allparamsPL,allcostsPL] = findPLparams(expNameunique,resultsfolder,length(paramNames),loadParameters,{paramname},paramNames);

    if ~foundPL
        disp('OBS no prev PL params found!!!')
        break
    end

    nanvals = isnan(allcostsPL);
    pvals = allparamsPL(~nanvals,strcmp(paramNames,paramname));
    costs = allcostsPL(~nanvals);
    [pvals,sind] = sort(pvals);
    costs = costs(sind,:);

    %remove double pvalues
    [pvalsu,ia,ic]  = unique(pvals);
    costsu = nan(size(pvalsu));
    for indpval = 1:length(pvalsu)
        costsu(indpval) = min(costs(ic == indpval));
    end
    pvals = pvalsu;
    costs = costsu;

    vals = [pvals,costs];
    if sum(~nanvals)>2
        [peakCostvals,peakPvals] = findpeaks(costs,pvals,'MinPeakProminence',0.3);
        [peakCostvals2,peakPvals2] = findpeaks(-costs,pvals,'MinPeakProminence',0.3);
        [peakPvals,uind] = unique([peakPvals;peakPvals2]);
        peakCostvals =[peakCostvals;-peakCostvals2];
        peakCostvals = peakCostvals(uind);
    end

    if sum(~nanvals)<=2 || isempty(peakPvals)
        disp('no peaks found!')
        % fprintf(fid,'no peaks found!\n');
        % fclose(fid);
        break
    else
        testrange = peakPvals';
        allcosts = vals;
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
    i = 1;
    for pvalue = testrange
        fprintf('-----Fixed %s to %0.4f-----------\n',paramname,pvalue)
        % fprintf(fid,'-----Fixed %s to %0.4f-----------\n',paramname,pvalue);
        allparams(indFixedParam) = pvalue;
        %check if a paramset already optimized, and use as start guess
        savefolder = sprintf('%s/%s_%d',PLfolder,paramname,pvalue);
        bestcostAll = peakCostvals(i);
        for rep = 1:nRepeats
            % load startguess from previous results
            try
                [optParam,optcost] = findBestParams(savefolder,0,nParams);
                disp('Loaded startguess')
                % fprintf(fid,'Loaded startguess');
            catch
                disp(['Couldnt load startguess from ' savefolder ', skipping this parameter value'])
                break
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
                allcosts((allcosts(:,1)==pvalue),2) = bestcostAll;
                allcosts((allcosts(:,1)==pvalue),1) = pvalue;
                save(sprintf('%s/allcostsLoaded_%s_%i_smoothed.mat',PLfolder,paramname,loopNumber),'allcosts') %save each round in case it craches or stops early
            end
        end

        %update the allcosts with the smoothed value
        allcosts((allcosts(:,1)==pvalue),2) = bestcostAll;
        allcosts((allcosts(:,1)==pvalue),1) = pvalue;

        save(sprintf('%s/allcostsLoaded_%s_%i_smoothed.mat',PLfolder,paramname,loopNumber),'allcosts') %save each round in case it craches or stops early

        i = i+1;

    end
    save(sprintf('%s/allcostsLoaded_%s_%i_smoothed.mat',PLfolder,paramname,loopNumber),'allcosts')

end


end

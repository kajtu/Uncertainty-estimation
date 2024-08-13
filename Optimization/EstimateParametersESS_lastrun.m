function [] =  EstimateParametersESS_lastrun(loopNumber,dateStart,experimentNames,experimentRange,numberOfSamples,startN)
% Same as EstimateParametersESS, but only one repetition and allowing for longer optimization time

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

%% Choose which subject to simulate based on the range of subjects and which loop number you are in
if ischar(experimentNames) || isstring(experimentNames)
    experimentNames = {experimentNames};
end

if strcmp(experimentNames{1}(1:3),'sim')% simulated data
    % Sample data
    n = mod(loopNumber,numberOfSamples)+1 + startN;
    expName = experimentNames{1};
    [data,~] = loadSampledMeasurementError(n,expName);
    data.(expName) = data; %create struct
    simulatedData = 1;
    disp(['-----------------------ESS lastrun. Exp: ' expName num2str(n) ' ---------------------------------'])
else % "measured" or other data
    simulatedData = 0;
    % Choose which experiment to simulate based on the range of experiments and which loop number you are in
    % Load the data to get all experiment names
    load('data.mat','data')
    % select experiments
    [~,experimentNames] = extractData(experimentNames,experimentRange,data);

    num = mod(loopNumber,length(experimentNames))+1;
    expName = experimentNames{num};
    disp(['----------------------- ESS lastrun. Exp: ' expName ' ---------------------------------'])
end



%% ESS options
% MEIGO OPTIONS I (COMMON TO ALL SOLVERS):
opts.ndiverse     = 'auto'; %100; %500; %5; %
if strcmp(data.meta.simulationFunction,'simulate_avatar_correction')
    opts.maxtime      = 5*3600; % MAX-Time of optmization, i.e how long the optimization will last %3600=1h %Maximum CPU time in seconds (Default 60)
else
    opts.maxtime      = 3600; % MAX-Time of optmization, i.e how long the optimization will last %3600=1h %Maximum CPU time in seconds (Default 60)
end
opts.maxeval      = 1e8; % max number of evals, i.e cost function calls
opts.log_var      = [];  %skip this

opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; % local solver, fmincon works in my experience best
opts.local.finish = opts.local.solver; %uses the local solver to check the best p-vector
opts.local.bestx = 0; % read the documentation, think it's best to leave at zero for now.
problem.f   = 'costFunction';%'costFunction_forchi2'; % %name of cost-function
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

%% Load parameters and data
resultsfolder = fullfile(basefolder,'Parameters');
[~,data,~,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues,~,~,~,~,~] = setup_simulations({expName},data,resultsfolder,1);

% settings for the cost function:
doPlot = 0;
dispErrors = 0;

%% Set bounds (separate for each subject, since paramdata is individual)
[lbOrig,ubOrig] = loadParamBounds(ind,data.parameters);
lbOrig = lbOrig(data.meta.paramsToOptimize);
ubOrig = ubOrig(data.meta.paramsToOptimize);
lb = log10(lbOrig);
ub = log10(ubOrig);
problem.x_L       = lb; % essOPT uses a problem structure where crucial information is specified
problem.x_U       = ub;

%% Find startguess
loadsummaryfile=0;
if simulatedData
    folders = dir(fullfile(resultsfolder,'ESS',['E_' expName num2str(n) '_*']));
else
    folders = dir(fullfile(resultsfolder,'ESS',['E_' expName '_*']));
end
[~,latestfolderInd] = max([folders.datenum]);
foldernames = {folders.name};
folderpaths = {folders.folder};
currentResultsFolder = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
fprintf('CurrentResultsFolder: %s\n',currentResultsFolder)

[optParam,loadedoptcost,~,loadedconstants] = findBestParams(currentResultsFolder,loadsummaryfile,length(paramNames)); %load from current results folder
fprintf('Estimateparameters: Loaded startguess from current results folder with cost %0.2f\n',loadedoptcost)
if sum(isnan(optParam))>0 % load from a previous optimization if current result is empty
    disp('Current results folder empty')
    [~,latestfolderInd] = max([folders.datenum]);
    latestfolderInd = latestfolderInd+1; %not the latest = this opt, but the one before
    if latestfolderInd > length(folders)
        latestfolderInd = latestfolderInd-2;
    end
    foldernames = {folders.name};
    folderpaths = {folders.folder};
    loadfolderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    [optParam,loadedoptcost,~,loadedconstants] = findBestParams(loadfolderName,loadsummaryfile,length(paramNames));
    fprintf('Estimateparameters: Loaded startguess with cost %0.2f from %s\n',loadedoptcost,loadfolderName)
end

if size(optParam) ~= size(lb)
    optParam = optParam';
end

optParam = max(optParam,lbOrig);%make sure the startguess is within the bounds
optParam = min(optParam,ubOrig);%make sure the startguess is within the bounds

optParam = log10(optParam);
optParam = max(optParam,lb);%make sure the startguess is within the bounds
optParam = min(optParam,ub);%make sure the startguess is within the bounds

problem.x_0=optParam;
allparams = origParamvalues;

%% Run optimization
format long
format compact
warning('off','all') % AMICI prints error-messages for all simulations (if they fail) so this will fix that annoying orange text form appearing
save(sprintf('%s/testsave-%i%s-%s.mat',currentResultsFolder,loopNumber,'lastrun',dateStart) ,'dateStart')

%% Solve with ESS
optim_algorithm = 'ess'; % 'multistart'; %  'cess'; % ESS IS BEST EVAH
Results = MEIGO(problem,opts,optim_algorithm,constants,allparams,simOptions,data,ind,dispErrors,doPlot); % Run the optimization
%% Save results
fitting_speed     = Results.time(end);
bestcost          = Results.fbest;
optParam          = Results.xbest;

%calculate cost to use for chi2test
costChi2 = costFunction_chi2(optParam,constants,allparams,simOptions,data,ind,dispErrors,doPlot);

% Save the results.
optParam = 10.^(optParam); % scale bestparam back to it's original size
bounds.lb = 10.^lb;
bounds.ub = 10.^ub;
dateEnd = datestr(now,'yymmdd-HHMMSS');
save(sprintf('%s/opt-ESS(%.3f)-%i%s.mat',currentResultsFolder,bestcost,loopNumber,'lastrun') ,'optParam','bestcost','costChi2','bounds','constants','ind','paramNames','simOptions','dateStart','dateEnd','fitting_speed','s','Results','problem')

warning ('on','all'); % yay error messages again


end

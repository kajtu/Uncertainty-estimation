function [] =  EstimateParametersESS_realdata(loopNumber,dateStart,dataName,randomstart)
% dataName: for example 'dataP1', name of the file containing the data for
% the wanted subject
% loopNumber: number in the loop calling this script (if not a loop, set to
% 1)
% dateStart: time and day when the script was called
% randomstart: use a random startguess or not

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

%% Load data
load(fullfile(basefolder,'Data',[dataName, '.mat']),'estimationData')
data = estimationData;
disp(dataName)

%%
if randomstart
    disp('Random start')
    randomstring = 'random';
    numRepeats = 5; % number of repetitions of running ESS optimizaton - don't run too much in case you are stuck in a terrible solution
else
    disp('No random start')
    randomstring = 'loaded';
    numRepeats = 20; % number of repetitions of running ESS optimizaton - run until not improving enough
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
problem.f   = 'costFunction_chi2';%name of cost-function
opts.iterprint = 1;

% MEIGO OPTIONS II (FOR ESS AND MULTISTART):
opts.local.iterprint = 1; % prints what going on during optimization

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
datain.(dataName) = data;
[~,data,~,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues,~,~,~,~,~] = setup_simulations({dataName},datain,resultsfolder,1);

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
folders = dir(fullfile(resultsfolder,'ESS',['E_' dataName '_*']));

if isempty(folders)
    mkdir(fullfile(resultsfolder,'ESS',['E_' dataName '_2024']))
    folders = dir(fullfile(resultsfolder,'ESS',['E_' dataName '_*']));
end

[~,latestfolderInd] = max([folders.datenum]);
foldernames = {folders.name};
folderpaths = {folders.folder};
currentResultsFolder = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
disp(['CurrentResultsFolder: ',currentResultsFolder])

if ~randomstart
    [optParam,loadedoptcost,~,loadedconstants] = findBestParams(currentResultsFolder,loadsummaryfile,length(paramNames)); %load from current results folder
end
if ~randomstart && sum(isnan(optParam))==0 %if there are current results and random start is not applied 
    fprintf('Estimateparameters: Loaded startguess from current results folder with cost %0.2f\n',loadedoptcost)
    optParam = optParam.*(rand(1, length(optParam))./20+0.975); %add small random change in the loaded parametes (+-2.5% at max)
else 
    if randomstart || isempty(folders) || length(folders) < 2 %load random start guess
        nParams = length(lb);
        optParam=(ubOrig-lbOrig).*rand(1, nParams)+lbOrig;
        disp('Estimateparameters: Loaded random startguess')
%         optParam = origParamvalues.*(rand(1, length(origParamvalues))./8+0.9375);
%         disp('Estimateparameters: Loaded random startguess based on origParamvalues')
    else % load from a previous optimization
        [~,latestfolderInd] = max([folders.datenum]);
        latestfolderInd = latestfolderInd+1; %not the latest = this opt, but the one before
        if latestfolderInd > length(folders)
            latestfolderInd = latestfolderInd-2;
        end
        foldernames = {folders.name};
        folderpaths = {folders.folder};
        loadfolderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
        [optParam,loadedoptcost,~,loadedconstants] = findBestParams(loadfolderName,loadsummaryfile,length(paramNames));
        fprintf('Estimateparameters: Loaded paramvalues with cost %0.2f from %s\n',loadedoptcost,loadfolderName)
        optParam = optParam.*(rand(1, length(optParam))./8+0.9375);% +-6.25% at max
    end
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
for r = 1:numRepeats
    fprintf(['----- ' dataName num2str(loopNumber) 'repetition num' num2str(r) ' of ' num2str(numRepeats) ' -----------------------\n'])
    
    %% Solve with ESS
    optim_algorithm = 'ess';
    Results = MEIGO(problem,opts,optim_algorithm,constants,allparams,simOptions,data,ind,dispErrors,doPlot); % Run the optimization
    %% Save results
    fitting_speed     = Results.time(end);
    bestcost          = Results.fbest;
    optParam          = Results.xbest;
    
    %update the start guess
    problem.x_0=optParam;

    %calculate cost to use for chi2test
    costChi2 = costFunction_chi2(optParam,constants,allparams,simOptions,data,ind,dispErrors,doPlot);
    
    % Save the results.
    optParam = 10.^(optParam); % scale bestparam back to it's original size
    bounds.lb = 10.^lb;
    bounds.ub = 10.^ub;
    dateEnd = datestr(now,'yymmdd-HHMMSS');
    save(sprintf('%s%sopt-ESS(%.3f)-%i%s-r%i.mat',currentResultsFolder,filesep,bestcost,loopNumber,randomstring,r) ,'optParam','bestcost','costChi2','allparams','bounds','constants','ind','paramNames','simOptions','dateStart','dateEnd','fitting_speed','s','Results','problem')
    
    %exit after 2 rounds if cost is not improving enough or after 1 round if much
    %worse than the loaded cost
    if ~randomstart && bestcost - loadedoptcost > 5
        break
    elseif r>1 && randomstart && bestcost-prevcost > -0.1
        break
    elseif r>1 && ~randomstart && (bestcost - loadedoptcost > 0.01 || bestcost-prevcost > -0.01)
        break
    else
        prevcost = bestcost;
    end
end
warning ('on','all'); % yay error messages again


end

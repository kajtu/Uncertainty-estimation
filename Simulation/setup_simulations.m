function [experimentNames,data,paramValues,constants,...
    paramNames,constantsNames,ynames,xnames,simulationoptions,inds,...
    origParamvalues,units,loadedcosts,meanParams,medianParams,paramuncertainty,allparams,okChi2] =setup_simulations(experimentNames,data,resultsFolder,doOptimization,experimentNameBase)


if ~doOptimization %if optimization, setup is run outside this function
    %make sure that the setup script can be found
    basefolder = split(pwd,'Uncertainty-estimation');
    addpath(fullfile(basefolder{1},'Uncertainty-estimation'));
    Setup('',0)
end
%names in the model
ynames = {'P_Aortic','Pperipheral','pressGrad MV','P Dmv','mv open','P Dav','av open','pressgrad AV','Ela','Elv','Qcaa','Qpc','P pulmvein','Qpvc','qLA','qLV','pLA','pLV','aaCorr','avCorr','mvCorr','pvCorr','Vla','Vlv'};
xnames = {'Ppvc','Qpv','Vla','Qmv','Vlv','Qav','Paa','Qaa','Ppc','mv_open','av_open'};

numExperiments = length(experimentNames);

%% Load parameters
nconstants = 17-3;
nparams = 28+4;
origParamvalues= nan(nparams,length(experimentNames));
constants = nan(nconstants,length(experimentNames));
inds = cell(length(experimentNames),1);
for e = 1:length(experimentNames)
    % Load constants and parameters
    [origParamvalues(:,e),constants(:,e),paramNames,constantsNames,units,inds{e}] = loadParameters(data.(experimentNames{e}));
end

if doOptimization
    paramValues = origParamvalues;
    loadedcosts = NaN;
    meanParams=NaN;
    paramuncertainty=NaN;
    medianParams = NaN;
    allparams =NaN;
else
    fprintf('SETUP: Loading all parameters...\n')
    if nargin > 4
        [experimentNames,paramValues,loadedcosts,meanParams,paramuncertainty.allokParams,...
            paramuncertainty.allokCosts,paramuncertainty.minValuesparams,...
            paramuncertainty.maxValuesparams,paramuncertainty.numincluded,...
            data,inds,constants,medianParams,allparams,okChi2] = findAllParams(resultsFolder,experimentNames,paramNames,data,inds,constants,experimentNameBase);
    else
        [experimentNames,paramValues,loadedcosts,meanParams,paramuncertainty.allokParams,...
            paramuncertainty.allokCosts,paramuncertainty.minValuesparams,...
            paramuncertainty.maxValuesparams,paramuncertainty.numincluded,...
            data,inds,constants,medianParams,allparams,okChi2] = findAllParams(resultsFolder,experimentNames,paramNames,data,inds,constants);
    end
    fprintf('SETUP: Loaded paramvalues for %d experiments\n',length(experimentNames))
end

%% Simulation settings
% sensi: sensitivity order - 0 = no calculation of sensitivities (default)
% maxsteps: maximum number of integration steps
simulationoptions = amioption('sensi',0,'maxsteps',1e4);
% set simulation tolerances
simulationoptions.atol = 1e-16;
simulationoptions.rtol = 1e-8;

% Calculate initial conditions for each subject based on the parameter values and data
ICs = zeros(length(experimentNames),length(xnames));
for e = 1:length(experimentNames)
    [ICs(e,:),~] = iccalcs(paramValues(:,e),constants(:,e),paramNames,constantsNames,data.(experimentNames{e}),inds{e});
    data.(experimentNames{e}).IC = ICs(e,:)';
end

%% Format the function output
if numExperiments == 1 && doOptimization
    experimentNames = experimentNames{1};
    inds= inds{1};
    data = data.(experimentNames);
end

end
%% Main script to reproduce the results in the manuscript
close all
clear

disp('*** Re-creating result figures. This might take several minutes. ***')

saveFigs = 1;

%% Set up: Load data, optimization results, and settings
% Add dependencies to matlab path
addpath(['.' filesep 'Simulation'])
addpath(['.' filesep 'Optimization'])
addpath(['.' filesep 'Data'])
addpath(['.' filesep 'Requirements'])

% Load the 100 simulated bootstrapped datasets
experimentName = 'simRandomSystematic';%simRandomSystematic or simRandomOnly
experimentNames = cell(1,100);
for n = 1:length(experimentNames)
    [datan,~] = loadSampledMeasurementError(n,experimentName);
    expName = [experimentName num2str(n)];
    data.(expName) = datan; %create struct
    experimentNames{n} = expName;
end

% Load true data
usetrueData = 1; 
load('dataSimulated.mat','dataSimulated','simulatedDataTrue')
trueData = dataSimulated.(experimentName);
trueData.allSimulations = simulatedDataTrue.allSimulations;
trueData.allParameters = simulatedDataTrue.allParameters;

% Load simulaiton settings and optimizaton results for all experiments
basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');
resultsFolder = fullfile(basefolder,'Parameters');
[experimentNames,data,bestparams,constants,...
    paramNames,constantsNames,ynames,xnames,simulationOptions,inds,...
    origParamvalues,units,loadedcosts,meanParams,medianParams,paramuncertainty,allparams] =setup_simulations(experimentNames,data,resultsFolder,0,experimentName);

% Create a folder to save the resulting figures in
thispath = split(pwd,'Uncertainty-estimation');
if saveFigs
    plotFolderName = fullfile(thispath{1},'Uncertainty-estimation',['Resultfigures_',experimentName,'_',datestr(now,'yymmdd_HHMM')]);
    mkdir(plotFolderName)
else
    plotFolderName = NaN;
end

%% Plot the figures

disp('Plotting data uncertainty (Fig 4)...')
cd Data
plotMeasurementErrors(experimentName,1,plotFolderName)
cd ../Plotting

disp('Plotting data uncertainty (Fig 5)...')
plotSampledData(experimentName,1,plotFolderName)

disp('Plotting method evaluation (Fig 6)...')
[allPLs,identifiabilityTable,identifiabilityTablePercent] = plotProfileLikelihood(experimentName,experimentNames,paramNames,data,inds,units,plotFolderName,saveFigs);

disp('Plotting the biomarker uncertainties for the subjects in the clinical study (Fig 7)...')
plotPL_realdata(plotFolderName)

disp('Plotting predictions (Fig S2)...')
plotPredictionUncertaintiesMain(data,trueData,experimentName, experimentNames,...
    bestparams,paramuncertainty,constants,paramNames,ynames,xnames,simulationOptions,inds,...
    saveFigs,plotFolderName)
cd ..

fprintf('***\nDONE. The figures and tables are saved in %s.***\n',plotFolderName)

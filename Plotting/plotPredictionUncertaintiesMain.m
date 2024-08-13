function [] = plotPredictionUncertaintiesMain(data,trueData,experimentName, experimentNames,...
    paramValues,paramuncertainty,constants,paramNames,ynames,xnames,simulationOptions,inds,...
    saveFigs,plotFolderName)

%% Set colors
l = length(experimentNames);
darkpurple = [0.9 0.5 0.9].*0.4;
lightpurple = [1 0.7 1];
purplegradient = [linspace(lightpurple(1),darkpurple(1),l)', linspace(lightpurple(2),darkpurple(2),l)', linspace(lightpurple(3),darkpurple(3),l)'];
colors = num2cell(purplegradient,2);

%% uncertainty simulations 
loadResults=1;
saveResults=1;
[minmaxSims] = findSimulationUncertainty(experimentNames,...
    paramuncertainty.allokParams,constants,inds,data,...
    simulationOptions,paramValues,saveResults,loadResults,experimentName);

%% Plot
plotPredictionUncertaintyCurves({minmaxSims},ynames,'FigS2_predictionUnc_curves',colors,trueData)

%% Save
if saveFigs
    savePDF = 1;
    saveAllFigures(plotFolderName,saveFigs,savePDF)
end
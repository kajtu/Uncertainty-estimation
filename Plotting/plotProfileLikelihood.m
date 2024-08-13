function [allPLs,identifiabilityTable,identifiabilityTablePercent] = plotProfileLikelihood(experimentName,experimentNames,paramNames,data,inds,units,plotFolderName,saveFigs)
% Plots the results from the profile likelihood for the simulated bootstrap
% data

%% Setup
basefolder = split(pwd,'4d-flow-mri-based-uncertainty');
basefolder = fullfile(basefolder{1},'4d-flow-mri-based-uncertainty');
resultsFolder = fullfile(basefolder,'Parameters');

if nargin <1
    close all
    saveFigs = 0;
    experimentName = 'simRandomOnly'; %simRandomSystematic or simRandomOnly

    addpath(genpath(basefolder))
    plotFolderName = fullfile(basefolder,'Plots','ParameterUncertainty',[experimentName datestr(now,'yymmdd-HHMM')]);

    % Load simulated data
    addpath '../Optimization'
    experimentNames = cell(1,100);
    for n = 1:length(experimentNames)
        [datan,~] = loadSampledMeasurementError(n,experimentName);
        expName = [experimentName num2str(n)];
        data.(expName) = datan; %create struct
        experimentNames{n} = expName;
    end

    addpath '../Simulation'
    [experimentNames,data,bestparams,constants,...
        paramNames,constantsNames,ynames,xnames,simulationOptions,inds,...
        origParamvalues,units,loadedcosts,meanParams,medianParams,paramuncertainty,allparams] =setup_simulations(experimentNames,data,resultsFolder,0,experimentName);
end

% Load simulated data
addpath '../Optimization'
experimentNames = cell(1,100);
for n = 1:length(experimentNames)
    [datan,errorn] = loadSampledMeasurementError(n,experimentName);
    expName = [experimentName num2str(n)];
    data.(expName) = datan; %create struct
    data.(expName).e = errorn;
    experimentNames{n} = expName;
end

% Load true data
load('dataSimulated.mat','dataSimulated','simulatedDataTrue')
trueData = dataSimulated.(experimentName);
trueData.allParameters = simulatedDataTrue.allParameters;


%% Load previously saved results
if strcmp(experimentName,'simRandomOnly')
    load('foundPL_simRandomOnly_240516-1421','foundPL','allPLs')
else
    load('foundPL_simRandomSystematic_240602-1053.mat','foundPL','allPLs')
end

%% OR re-create them
% thresholdChi2all = zeros(size(experimentNames));
% allPLs = cell(size(experimentNames));
% foundPL = logical(size(experimentNames));
% % load all PLs
% for e = 1:length(experimentNames)
%     experiment = experimentNames{e};
%     [dgf_data,~] = degreesOfFreedom(data.(experimentNames{e}),inds{e});
%     thresholdChi2all(e) = chi2inv(0.95,dgf_data);
%
%     % load all PLs for this experiment
%     loadParameters=0;
%     [foundPL(e),allPLs{e}] = findPLparams(experiment,resultsFolder,length(paramNames),loadParameters,paramNames);
% end
%
% savename = sprintf('foundPL_%s_%s.mat',experimentName,datestr(now,'yymmdd-hhMM'));
% save(savename,'foundPL','allPLs','thresholdChi2all','dgf_data')
% disp(['Saved last PL results in: ' savename])

%% correct params that could not be optimized to its simulation value
for e = 1:length(experimentNames)
    % offset correction parameters
    %-40 since the bounds are -40 to 40 but opt bounds are 0 to 80 = -40+40 to 40+40.
    existingparams = fieldnames(allPLs{e});
    if sum(strcmp(existingparams,'mvCorr'))>0
        allPLs{e}.('mvCorr')(:,1) = allPLs{e}.('mvCorr')(:,1) - 40;
    end
    if sum(strcmp(existingparams,'onset_LV'))>0
        %onset_LV: range set to 1-2 instead of -0.5 to 0.5 --> take onset_LV-1.5
        allPLs{e}.('onset_LV')(:,1) = allPLs{e}.('onset_LV')(:,1) - 1.5;
    end
end

%same for the true values:
% fix the onset LV
pind = ismember(paramNames,'onset_LV');
trueData.allParameters(pind) = trueData.allParameters(pind) - 1.5;

%% plot all PLs
l = length(experimentNames);
purple = [0.9 0.5 0.9];
darkpurple = [0.9 0.5 0.9].*0.4;
lightpurple = [1 0.7 1];
purplegradient = [linspace(lightpurple(1),darkpurple(1),l)', linspace(lightpurple(2),darkpurple(2),l)', linspace(lightpurple(3),darkpurple(3),l)'];
green = [0 0.55 0.35];

%% physiological order of parameters
PVparams = {'Ppu','Rpu','Cpvc','Lpv','Rpv'};
LAparams = {'Emax_LA','Emin_LA','k_diast_LA','k_syst_LA','m1_LA','m2_LA','onset_LA'};
MVparams = {'Lmv','Rmv'};
LVparams= {'Emax_LV','Emin_LV','k_diast_LV','k_syst_LV','m1_LV','m2_LV','onset_LV'};
AVparams = {'Lav','ELCo'};
Aortaparams = {'Caa','Lao','Rao'};
Perparams = {'Rtot','Ctot'};
paramNamesPhysOrder = [PVparams,LAparams,MVparams,LVparams,AVparams,Aortaparams,Perparams];
paramNamesPhysOrder = [paramNamesPhysOrder,{'mvCorr','avCorr','aaCorr','pvCorr'}];

%% setup variables
if strcmp(experimentName,'simRandomSystematic')
    paramNamesPlot =  fixUnderscore(paramNamesPhysOrder(1:end-4));
else
    paramNamesPlot =  fixUnderscore(paramNamesPhysOrder);
end

identifiable = zeros(length(paramNamesPlot),length(experimentNames));
onlylowerBoundFound = zeros(length(paramNamesPlot),length(experimentNames));
onlyupperBoundFound = zeros(length(paramNamesPlot),length(experimentNames));
unidentifiable = zeros(length(paramNamesPlot),length(experimentNames));

ub = zeros(length(paramNamesPlot),length(experimentNames));
lb = zeros(length(paramNamesPlot),length(experimentNames));
foundTrue = zeros(length(paramNamesPlot),length(experimentNames));

%% Figure with all PLs
figure('Name',['FigS1_ProfileLikelihood_all_' experimentName])
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 26;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout('flow','TileSpacing','loose','Padding','compact')
% one plot for each param
for i = 1:length(paramNamesPlot)
    p = find(ismember(paramNames,paramNamesPhysOrder{i}));
    nexttile
    hold on
    xlabel([paramNamesPlot{i} ' (' units.param{p} ')'])
    ylabel('Obj. func. value')
    if ~strcmp(paramNames{p},'mvCorr')
        xline(trueData.allParameters(p),'--','color',green,'LineWidth',1.5);
    end
    for e = 1:length(experimentNames)
        pname = paramNames{p};
        if foundPL(e) &&  ismember(pname,fieldnames(allPLs{e}))
            %sort
            [~,sind]=sort(allPLs{e}.(pname)(:,1));

            %check identifiability
            costs = allPLs{e}.(pname)(sind,2);
            pvals = allPLs{e}.(pname)(sind,1);
            [mincost, midind] = min(costs);
            limit = mincost+chi2inv(0.95,1);
            aboveinds = find(costs > limit);%thresholdChi2

            %check identifiability
            if isempty(aboveinds) % if no crossings found
                identifiable(p,e) = 0;
                onlylowerBoundFound(p,e) = 0;
                onlyupperBoundFound(p,e) = 0;
                unidentifiable(p,e) = 1;
            elseif sum(aboveinds>midind)>0 && sum(aboveinds<midind)>0 %if found crossing right and left of midpoint
                identifiable(p,e) = 1;
                onlylowerBoundFound(p,e) = 0;
                onlyupperBoundFound(p,e) = 0;
                unidentifiable(p,e) = 0;
            elseif sum(aboveinds>midind)>0 %if found crossing right of midpoint
                identifiable(p,e) = 0;
                onlylowerBoundFound(p,e) = 0;
                onlyupperBoundFound(p,e) = 1;
                unidentifiable(p,e) = 0;
            elseif sum(aboveinds<midind)>0  %if found crossing left of midpoint
                identifiable(p,e) = 0;
                onlylowerBoundFound(p,e) = 1;
                onlyupperBoundFound(p,e) = 0;
                unidentifiable(p,e) = 0;
            end

            % find bounds
            ub(p,e) = max(pvals(costs <= limit));
            lb(p,e) = min(pvals(costs <= limit));

            %if mvccorr, the true value is different for each sampled data
            if strcmp(pname,'mvCorr')
                trueparamvalue = data.(experimentNames{e}).e.mv.systematic;
            else
                trueparamvalue =trueData.allParameters(p);
            end

            %check if found true value
            if trueparamvalue<= ub(p,e) && trueparamvalue >= lb(p,e)
                foundTrue(p,e) = 1;
                plot(allPLs{e}.(pname)(sind,1),allPLs{e}.(pname)(sind,2),'*-','color',purplegradient(e,:),'LineWidth',1,'MarkerSize',3)
            else
                foundTrue(p,e) = 0;
                plot(allPLs{e}.(pname)(sind,1),allPLs{e}.(pname)(sind,2),'k*-','LineWidth',1,'MarkerSize',3)
            end

            if max(allPLs{e}.(pname)(sind,2)) > 500
                ylim([0 500])
            end

        end

    end
end

%% Summary figure for publication
l = 100;
white = [1 1 1];
green = [0 0.55 0.35];
greengradient = [linspace(white(1),green(1),l)', linspace(white(2),green(2),l)', linspace(white(3),green(3),l)'];

purple = [0.9 0.5 0.9];
darkpurple = [0.9 0.5 0.9].*0.4;
lightpurple = [1 0.7 1];
purplegradient = [linspace(lightpurple(1),darkpurple(1),l)', linspace(lightpurple(2),darkpurple(2),l)', linspace(lightpurple(3),darkpurple(3),l)'];

selectedParamNames = {'k_diast_LA','m2_LV','Caa'};
selectedParamNamesPlot = {'k_d_i_a_s_t__L_A','m2_L_V','Caa'};
figure('Name',['Fig6_methodevaluation_summary_' experimentName])
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 15;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(3,3,'TileSpacing','loose','Padding','compact')

foundTrueParams = sum(foundTrue,2);
[~,paramorder] = ismember(paramNamesPhysOrder,paramNames);
if strcmp(experimentName,'simRandomSystematic')
    [~,paramorder] = ismember(paramNamesPhysOrder(1:end-4),paramNames(1:end-4));
else
    selectedParamNames = {'m2_LV','Caa'};%{'Caa','Rao','m2_LV'};
    selectedParamNamesPlot = {'m2_L_V','Caa'};%{'Caa','Rao','m2_L_V'};
end

nexttile([3,1])
hHM = heatmap({'% of samples'},paramNamesPlot,100*foundTrueParams(paramorder)./length(experimentNames),'colormap',greengradient,'ColorLimits',[0 100]);
set(gca,'FontSize',10)
title('True value within conf. interval?')

for i = 1:length(selectedParamNamesPlot)
    % bounds
    nexttile
    if i==1
        title('Conf. interval')
    end
    p = find(strcmp(paramNames,selectedParamNames{i}));
    hold on
    for e = 1:length(experimentNames)
        if foundTrue(p,e)
            plot([lb(p,e),ub(p,e)],[e,e],'.-','linewidth',1.5,'color',purplegradient(e,:),'MarkerSize',5)
        else
            plot([lb(p,e),ub(p,e)],[e,e],'k.-','linewidth',1.5,'MarkerSize',5)
        end
    end
    xline(trueData.allParameters(p),'--','LineWidth',1.2,'color',[0 0.5 0]);
    xlabel([paramNamesPlot{ismember(paramNamesPlot,paramNames(p))} ' (value)'])
    ylabel('Sampled data n.')

    % full PL
    nexttile
    if i==1
        title('Full profile')
    end
    hold on
    for e = 1:length(experimentNames)
        pname = selectedParamNames{i};
        %sort
        [~,sind]=sort(allPLs{e}.(pname)(:,1));
        %check if found true value
        if foundTrue(p,e)
            okparamp = plot(allPLs{e}.(pname)(sind,1),allPLs{e}.(pname)(sind,2),'.-','color',purplegradient(e,:),'LineWidth',1,'MarkerSize',7);
        else
            notokparamp = plot(allPLs{e}.(pname)(sind,1),allPLs{e}.(pname)(sind,2),'k.-','LineWidth',1,'MarkerSize',7);
        end
    end
    truep = xline(trueData.allParameters(p),'--','LineWidth',1.2,'color',[0 0.5 0]);
    xlabel([paramNamesPlot{ismember(paramNamesPlot,paramNames(p))} ' (value)'])
    ylabel('Obj. func value')
end
legend([okparamp,notokparamp,truep],sprintf('Found true\nvalue'),sprintf('Didn''t find\n true value'),'True value',...
    'Position',[0.822 0.173 0.149 0.14])

%% Create tables
numExperimentsFound = sum([sum(identifiable,2),sum(onlylowerBoundFound,2),sum(onlyupperBoundFound,2),sum(unidentifiable,2)],2);
indsinclude = numExperimentsFound>0;%only include params where we found PL

if strcmp(experimentName,'simRandomSystematic')
    paramnames = paramNames(1:end-4);
    trueparams = trueData.allParameters(1:end-4)';
    [~,paramorder] = ismember(paramNamesPhysOrder(1:end-4),paramnames);
else
    paramnames = paramNames;
    trueparams = trueData.allParameters';
    paramorder = 1:sum(indsinclude);
end
paramnamesPlot = paramnames(paramorder);


lbtab = lb(paramorder,:);
ubtab = ub(paramorder,:);
trueparams = trueparams(paramorder);

lowertext = split(sprintf('%0.3f +- %0.3f (%0.3f)x',[mean(lbtab,2),std(lbtab,[],2),min(lbtab,[],2)]'),'x');
uppertext = split(sprintf('%0.3f +- %0.3f (%0.3f)x',[mean(ubtab,2),std(ubtab,[],2),max(ubtab,[],2)]'),'x');
trueparams = split(sprintf('%0.3fx',trueparams),'x')
boundsTable = table(lowertext(1:end-1),uppertext(1:end-1),trueparams(1:end-1),...
    'Rownames',paramnamesPlot,'VariableNames',{'Lower bound: mean +-sd (min)','Upper bound: mean +-sd (max)','True value'})

identifiabilityTable = table(sum(foundTrue(indsinclude,:),2), ...
    sum(identifiable(indsinclude,:),2),...
    sum(onlylowerBoundFound(indsinclude,:),2),...
    sum(onlyupperBoundFound(indsinclude,:),2),...
    sum(unidentifiable(indsinclude,:),2),...
    'Rownames',paramnamesPlot,...
    'VariableNames',{'Found true','Identifiable','Semi-identifiable: only lower bound found','Semi-identifiable: only upper bound found','Not identifiable'})

n  =numExperimentsFound(indsinclude);
identifiabilityTablePercent = table(100*sum(foundTrue(indsinclude,:),2)./n,... ...
    100*sum(identifiable(indsinclude,:),2)./n,...
    100*sum(onlylowerBoundFound(indsinclude,:),2)./n,...
    100*sum(onlyupperBoundFound(indsinclude,:),2)./n,...
    100*sum(unidentifiable(indsinclude,:),2)./n,...
    'Rownames',paramnamesPlot,...
    'VariableNames',{'Found true(%)','Identifiable (%)','Semi-identifiable: only lower bound found (%)','Semi-identifiable: only upper bound found (%)','Not identifiable (%)'})

%sort
identifiabilityTable = identifiabilityTable(paramorder,:)
identifiabilityTablePercent = identifiabilityTablePercent(paramorder,:)
boundsTable = boundsTable(paramorder,:)


%% Save and close all figures
if saveFigs
    if nargin <1
        mkdir(plotFolderName)
    end
    savePDF = 1;
    saveAllFigures(plotFolderName,saveFigs,savePDF)
    writetable(identifiabilityTablePercent,fullfile(plotFolderName,'identifiabilityTablePercent.xlsx'),"WriteRowNames",1)
    writetable(identifiabilityTable,fullfile(plotFolderName,'identifiabilityTable.xlsx'),"WriteRowNames",1)
    writetable(boundsTable,fullfile(plotFolderName,'boundsTable.xlsx'),"WriteRowNames",1)
end


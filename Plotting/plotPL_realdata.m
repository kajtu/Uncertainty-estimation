function plotPL_realdata(plotFolderName)
% Script to plot results from the parameter estimation method applied to a
% clinical dataset from the Health study (Tunedal 2023).

saveFigs = 1;
basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');

if nargin <1
    close all
    date = datestr(now, 'yymmdd');
    plotFolderName = fullfile(basefolder,'Plots',['clinicalexample_' date]);
    mkdir(plotFolderName)

    % add all project folders to the matlab path
    addpath(genpath(fullfile(basefolder,'Optimization')))
    addpath(genpath(fullfile(basefolder,'Data')))
    addpath(genpath(fullfile(basefolder,'Modelfiles')))
    addpath(genpath(fullfile(basefolder,'Requirements')))
    addpath(genpath(fullfile(basefolder,'Simulation')))

    % setup AMICI and MEIGO toolboxes
    run(fullfile(basefolder, 'Requirements', 'AMICI-0.10.11_SS_eventFix', 'matlab', 'installAMICI.m'))
    run(fullfile(basefolder, 'Requirements', 'MEIGO', 'install_MEIGO.m'))

    loadresults = 1;
else
    loadresults = 0;
end

%% Setup which biomarkers, units, and which subjects to plot
paramNamesPlot = {'Caa', 'Rao','m2_LV','k_syst_LV','Emax_LA','Cpvc'};
paramNamesPlotNice = {'Caa', 'Rao','m2_L_V','ksystLV','EmaxLA','Cpvc'};


experimentNames = {'E_dataP78','E_dataP1','E_dataP3','E_dataP24','E_dataP36','E_dataP33'};
plotNames = {'Control 1','Control 2','Control 3','T2D+HT 1','T2D+HT 2','T2D+HT 3'};

units.param = {'ml/mmHg','mmHg*s/ml','mmHg*s/ml','mmHg*s^2/ml','mmHg*s/ml',...
    'ml/mmHg','cm^2','ml/mmHg','ml/mmHg','ml/mmHg','ml/mmHg','ml/mmHg',...
    'mmHg*s^2/ml','mmHg*s^2/ml','mmHg*s^2/ml','mmHg','mmHg*s/ml','mmHg*s/ml',...
    '-','-','-','-','-','-','-','-','s','s','ml','ml','ml','ml'};

paramNames = {'Cpvc' 'Rpu' 'Rpv' 'Lpv' 'Rtot' 'Ctot' 'ELCo' 'Caa' 'Emax_LA' 'Emax_LV' 'Emin_LA' 'Emin_LV' 'Lao' 'Lav'...
    'Lmv' 'Ppu' 'Rao' 'Rmv' 'k_diast_LA' 'k_diast_LV' 'k_syst_LA'...
    'k_syst_LV' 'm1_LA' 'm1_LV' 'm2_LA' 'm2_LV' 'onset_LA' 'onset_LV'};


%% find PL results for all datasets
if loadresults
    %% Load all parameters
    resultsFolder = fullfile(basefolder,'Parameters');
    allPLs = cell(size(experimentNames));
    inds = cell(size(experimentNames));
    simOptions = cell(size(experimentNames));
    foundPL = logical(size(experimentNames));
    bestparam = nan(length(experimentNames),length(paramNames));
    bestcost = nan(length(experimentNames),1);
    thresholdChi2all = nan(length(experimentNames),1);
    % load all PLs
    for e = 1:length(experimentNames)
        experiment = experimentNames{e}(3:end);
        experimentNames{e} = experiment;
        disp(experiment)

        %load data
        load(fullfile(basefolder,'Data',[experiment, '.mat']),'estimationData')
        data.(experiment) = estimationData;

        % setup
        resultsfolder = fullfile(basefolder,'Parameters');
        [~,dataforcost{e},~,constants{e},~,constantsNames,ynames,xnames,simOptions{e},inds{e},origParamvalues{e},~,~,~,~,~] = setup_simulations({experiment},data,resultsfolder,1);

        [dgf_data,~] = degreesOfFreedom(data.(experiment),inds{e});
        thresholdChi2all(e) = chi2inv(0.95,dgf_data);

        % load all PLs for this experiment
        loadParameters=1;
        [foundPL(e),allPLs{e},~,~,bestparam(e,:),bestcost(e)] = findPLparams(experiment,resultsFolder,length(paramNames),loadParameters,paramNames);
    end

    %% Find parameter bounds
    ub = zeros(length(paramNamesPlot),length(experimentNames));
    lb = zeros(length(paramNamesPlot),length(experimentNames));

    % one plot for each param
    for i = 1:length(paramNamesPlot)
        p = find(ismember(paramNames,paramNamesPlot{i}));
        for e = 1:length(experimentNames)
            pname = paramNames{p};
            if foundPL(e) &&  ismember(pname,fieldnames(allPLs{e}))

                %sort
                [~,sind]=sort(allPLs{e}.(pname)(:,1));
                costs = allPLs{e}.(pname)(sind,2);
                pvals = allPLs{e}.(pname)(sind,1);
                [mincost, midind] = min(costs);
                limit = mincost+chi2inv(0.95,1);

                % find bounds
                ub(p,e) = max(pvals(costs <= limit));
                lb(p,e) = min(pvals(costs <= limit));
            end

        end
    end
    save('loadedPLs_clinicaldataset.mat','lb','ub','bestparam','bestcost','thresholdChi2all')

else
    %% OR load them from previous results:
    load('loadedPLs_clinicaldataset.mat','lb','ub','bestparam','bestcost','thresholdChi2all')
end

%% Check fit to data
table(bestcost,thresholdChi2all,bestcost<thresholdChi2all, ...
    'Rownames',experimentNames,'VariableNames',{'Best cost','Threshold','Below threshold?'})


%% Calculate individual and cohort variability
percAll = nan(size(paramNamesPlot));
perc = nan(size(paramNamesPlot));
for i = 1:length(paramNamesPlot)
    p = find(ismember(paramNames,paramNamesPlot{i}));
    meanbestparam = mean(bestparam(:,p));
    middlep = lb(p,:)+ ((ub(p,:)-lb(p,:))./2);%mid point of the confidence interval
    sd = (ub(p,:)-middlep) / 1.96; %95% conf interval corresponds to 1.96 sd
    perc(i) = mean(100* (sd./meanbestparam ));
    sdbestparam = std(bestparam(:,p));
    percAll(i) = 100* (sdbestparam / meanbestparam );
end
[vals,sortinds] = sort(perc./percAll);
paramNamesPlot = paramNamesPlot(sortinds);
table(perc(sortinds)',percAll(sortinds)',perc(sortinds)'./percAll(sortinds)','rownames',paramNamesPlot,'variablenames',{'Individual variation','Cohort variation','Ratio'})

paramNamesPlotNice = paramNamesPlotNice(sortinds);
perc = perc(sortinds);
percAll = percAll(sortinds);

%% Check significance (overlapping 95% confidence intervals)
% as an indication of significant differences between subjects, we check if
% the 95% confidence intervals overlap or not
sign = cell(length(paramNamesPlot),1);
comps = {};
for i = 1:length(paramNamesPlot)
    signp = cell(length(experimentNames),length(experimentNames));
    p = find(ismember(paramNames,paramNamesPlot{i}));
    for e = 1:length(experimentNames)
        % compare to all other subjects
        allother = 1:length(experimentNames);
        allother(allother<=e) = [];%already tested
        middlep = lb(p,e)+ ((ub(p,e)-lb(p,e))/2);
        for c = 1:length(allother)
            comp = allother(c);
            if i == 1
                comps = [comps,{[e,comp]}];
            end
            if lb(p,e) > ub(p,comp) || ub(p,e) < lb(p,comp) %no overlapping confidence intervals, definetely different between these subjects
                signp{e,comp} = '**';
                % sign(comp,e) = '**';
            elseif middlep > ub(p,comp) || middlep < lb(p,comp) %less than half overlapping, probably different (~ p=0.05)
                signp{e,comp} = '*';
            else % overlapping more than half, probably not different ( ~ p > 0.05)
                signp{e,comp} = '-';
            end
        end
    end
    sign{i} = signp;
end

%% Plot results figure for publication
darkpurple = [0.9 0.5 0.9].*0.4;
lightpurple = [1 0.7 1];
l = length(experimentNames);
purplegradient = [linspace(lightpurple(1),darkpurple(1),l)', linspace(lightpurple(2),darkpurple(2),l)', linspace(lightpurple(3),darkpurple(3),l)'];

letters = 'A':'Z';
figure('Name','Fig7_ProfileLikelihood_clinicalexample')
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 12+5+3;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(4,6,'TileSpacing','loose','Padding','compact')

for i = 1:length(paramNamesPlot)
    p = find(ismember(paramNames,paramNamesPlot{i}));
    ax1=nexttile([1,3]);
    hold on
    ylabel([paramNamesPlotNice{i} ' (' units.param{p} ')'])    
    for e = 1:length(experimentNames)
        middlep = lb(p,e)+ ((ub(p,e)-lb(p,e))/2);
        sd = (ub(p,e)-middlep) / 1.96; %95% conf interval corresponds to 1.96 sd
        if e > 3 % patient
            errorbar(e,middlep,sd,'.','color',purplegradient(end,:),'linewidth',2)
            plot(e,bestparam(e,p),'o','color',purplegradient(end,:),'MarkerFaceColor',purplegradient(end,:),'MarkerSize',3)
        else % control
            errorbar(e,middlep,sd,'.','color',purplegradient(1,:),'linewidth',2)
            plot(e,bestparam(e,p),'o','color',purplegradient(1,:),'MarkerFaceColor',purplegradient(1,:),'MarkerSize',3)
        end
    end
    xticks(1:length(plotNames))
    xticklabels(plotNames)
    xlim([0.5,length(plotNames)+0.5])
    if min(lb(p,:))>10
        y1 = round(min(lb(p,:)),0);
    else
        y1 = round(min(lb(p,:)),3);
    end
    if max(ub(p,:))>10
        y2 = round(max(ub(p,:)),0);
    else
        y2 = round(max(ub(p,:)),3);
    end
    if y1 < min(lb(p,:))*0.95
        y1 = min(lb(p,:))*0.95;
    end
    yticks([y1 y2]);
    yd = 0.05*(length(comps)-1) + 1;
    ymax = max(ub(p,:))*(yd+0.02)*1.01;
    ymax = max(max(ub(p,:))*1.05,ymax);
    ylim([min(lb(p,:))*0.95,ymax])

    %significance
    for c = 1:length(comps)
        yd = 0.05*(c-1) + 1;
        if strcmp(sign{i}(comps{c}(1),comps{c}(2)),'**')
            plot(comps{c},[max(ub(p,:))*yd,max(ub(p,:))*yd],'k-')
            plot(mean(comps{c}),max(ub(p,:))*(yd+0.02),'k*','markersize',4,'linewidth',0.5)
        end
    end
    set(gca,'FontSize',9,'FontName','Calibri')
    title(letters(i),'FontSize',11,'FontName','Calibri')
    ax1.TitleHorizontalAlignment = 'left';

    if i ==1
        legend({'Standard deivation (sd)','Best fit to data','No overlapping CI'},'Location','Northwest')
    end
end
lettersend =i;

% percent plot
green = [0 0.55 0.35];
for i = 1:length(paramNamesPlot)
    p = find(ismember(paramNames,paramNamesPlot{i}));
    ax1=nexttile;
    hold on
    meanbestparam = mean(bestparam(:,p));
    middlep = lb(p,:)+ ((ub(p,:)-lb(p,:))./2);
    sd = (ub(p,:)-middlep) / 1.96; %95% conf interval corresponds to 1.96 sd
    percP = 100* (sd./meanbestparam );

    bar([perc(i),percAll(i)],'FaceColor',green)
    errorbar([perc(i),percAll(i)],[std(percP),0],'k.')

    ylabel(['% of mean ' paramNamesPlotNice{i}])
    xticks([1,2])
    xticklabels({'Individual sd','Cohort sd'})
    set(gca,'FontSize',9,'FontName','Calibri')
    title(letters(i+lettersend),'FontSize',11,'FontName','Calibri')
    ax1.TitleHorizontalAlignment = 'left';
end


%% Save and close all figures
if saveFigs
    savePDF = 1;
    saveAllFigures(plotFolderName,saveFigs,savePDF)
end
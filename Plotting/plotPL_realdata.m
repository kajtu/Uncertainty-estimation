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
plotNames = {'C1','C2','C3','P1','P2','P3'};

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
linespace=0.05;
darkpurple = [0.9 0.5 0.9].*0.4;
lightpurple = [1 0.7 1];
l = length(experimentNames);
purplegradient = [linspace(lightpurple(1),darkpurple(1),l)', linspace(lightpurple(2),darkpurple(2),l)', linspace(lightpurple(3),darkpurple(3),l)'];

letters = 'A':'Z';
figure('Name','Fig7_ProfileLikelihood_clinicalexample')
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 12+3; 
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(3,12,'TileSpacing','loose','Padding','compact')

%example
nexttile([1,3])
axis('off')

% results
for i = 1:length(paramNamesPlot)
    p = find(ismember(paramNames,paramNamesPlot{i}));
    if i == 4
        nexttile([1,3]);
        axis('off')
    end
    ax1=nexttile([1,3]);
    hold on
    ylabel([paramNamesPlotNice{i} ' (' units.param{p} ')'])    
    for e = 1:length(experimentNames)
        xl=xline(e,':','Color',[0.1 0.1 0.1],'linewidth',0.8);
        middlep = lb(p,e)+ ((ub(p,e)-lb(p,e))/2);
        sd = (ub(p,e)-middlep) / 1.96; %95% conf interval corresponds to 1.96 sd
        if e > 3 % patient
            pe=errorbar(e,middlep,sd,'.','color',darkpurple,'linewidth',1.1,'markersize',1);
            pb=plot(e,bestparam(e,p),'o','color',darkpurple,'MarkerFaceColor',darkpurple,'MarkerSize',3);
        else % control
            ce=errorbar(e,middlep,sd,'.','color',lightpurple.*0.9,'linewidth',1.1,'markersize',1);
            cb=plot(e,bestparam(e,p),'o','color',lightpurple.*0.9,'MarkerFaceColor',lightpurple.*0.9,'MarkerSize',3);
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
    s=strcmp(sign{i},'**');
    numsign = sum(s(:));
    yd = linespace*(numsign) + 1 + (linespace*0.9)*length(experimentNames);
    ymax = max(ub(p,:))*(yd+0.02)*1.01;
    ymax = max(max(ub(p,:))*1.05,ymax);
    ylim([min(lb(p,:))*0.95,ymax])

    %significance
    yd=1;
    startingc = 1;
    for c = 1:length(comps)
        if strcmp(sign{i}(comps{c}(1),comps{c}(2)),'**')
            if startingc ~= comps{c}(1)
                yd = yd+linespace;
            end
            startingc = comps{c}(1);
            yd = yd+linespace;
            psign=plot(comps{c},[max(ub(p,:))*yd,max(ub(p,:))*yd],'k-');
        end
    end
    set(gca,'FontSize',9,'FontName','Calibri')
    title(letters(i),'FontSize',11,'FontName','Calibri')
    ax1.TitleHorizontalAlignment = 'left';

    if i ==1
        legend([xl,cb,ce,pb,pe,psign],...
            {sprintf('C: control\nP: HT+T2D'),'Standard deivation (C)','Best fit to data (C)','Standard deivation (P)','Best fit to data (P)','No overlapping CI'},...
            'Position',[0.02 0.415931811460952 0.2 0.19],'box','off')
    end
end
lettersend =i;


% percent plot
green = [0 0.55 0.35];
for i = 1:length(paramNamesPlot)
    p = find(ismember(paramNames,paramNamesPlot{i}));
    ax1=nexttile([1 2]);
    hold on
    meanbestparam = mean(bestparam(:,p));
    middlep = lb(p,:)+ ((ub(p,:)-lb(p,:))./2);
    sd = (ub(p,:)-middlep) / 1.96; %95% conf interval corresponds to 1.96 sd
    percP = 100* (sd./meanbestparam );

    b1=bar(1,perc(i),'FaceColor',[184 84 184]./255);
    b2=bar(2,percAll(i),'FaceColor',green);

    e1=errorbar([perc(i),percAll(i)],[std(percP),0],'k.');

    ylabel(['% of mean ' paramNamesPlotNice{i}])
    xticks([])
    set(gca,'FontSize',9,'FontName','Calibri')
    title(letters(i+lettersend),'FontSize',11,'FontName','Calibri')
    ax1.TitleHorizontalAlignment = 'left';
end
legend([b1,e1,b2],{sprintf('Individual\nstandard\ndeviation\n '),sprintf('Variation\namong\nindividuals\n '),sprintf('Cohort\nstandard\ndeviation\n ')},'Position',[0.925 0.12 0.05 0.23],'box','off')

%% Save and close all figures
if saveFigs
    savePDF = 1;
    saveAllFigures(plotFolderName,saveFigs,savePDF)
end
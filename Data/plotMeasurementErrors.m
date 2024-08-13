function plotMeasurementErrors(experimentName,saveStuff,plotFolderName)
% create measurement error figure

if nargin <1
    close all
    clear
    saveStuff = 1;
    experimentName = 'simRandomSystematic'; %simRandomSystematic or simRandomOnly
end

% Add dependencies to matlab path
addpath(['..' filesep 'Requirements'])

% Load true data
load('dataSimulated.mat','dataSimulated','simulatedDataTrue')
trueData = dataSimulated.(experimentName);
trueParams = trueData.parameters;


%% Load simulated data
addpath '../Optimization'
experimentNames = cell(1,100);
Ctoterr = nan(length(experimentNames),1);
Rtoterr = nan(length(experimentNames),1);
Emax_LVerr = nan(length(experimentNames),1);

SV_MV = nan(length(experimentNames),1);
SV_AV = nan(length(experimentNames),1);
SV_AA = nan(length(experimentNames),1);
SV_PV = nan(length(experimentNames),1);

SV_MV_randonly = nan(length(experimentNames),1);
SV_AV_randonly = nan(length(experimentNames),1);
SV_AA_randonly = nan(length(experimentNames),1);
SV_PV_randonly = nan(length(experimentNames),1);

SV_MV_systonly = nan(length(experimentNames),1);
SV_AV_systonly = nan(length(experimentNames),1);
SV_AA_systonly = nan(length(experimentNames),1);
SV_PV_systonly = nan(length(experimentNames),1);

for n = 1:length(experimentNames)
    [datan,e] = loadSampledMeasurementError(n,experimentName);
    [Ctot,Rtot,Emax_LV] = calculateDatabasedParameters(datan);
    Ctoterr(n) = Ctot-trueData.parameters.Ctot.mean;
    Rtoterr(n) = Rtot-trueData.parameters.Rtot.mean;
    Emax_LVerr(n) = Emax_LV-trueData.parameters.Emax_LV.mean;

    SV_MV(n) = trapz(datan.MV.time,datan.MV.mean);
    SV_AV(n) = trapz(datan.AV.time,datan.AV.mean);
    SV_AA(n) = trapz(datan.AA.time,datan.AA.mean);
    SV_PV(n) = trapz(datan.PV.time,datan.PV.mean);


    %add only random errors and calculate stroke volume
    SV_MV_randonly(n) =  trapz(datan.MV.time, trueData.MV.mean + e.mv.random);
    SV_AV_randonly(n) =  trapz(datan.MV.time, trueData.AV.mean + e.av.random);
    SV_AA_randonly(n) =  trapz(datan.MV.time, trueData.AA.mean + e.aa.random);
    SV_PV_randonly(n) =  trapz(datan.MV.time, trueData.PV.mean + e.pv.random);

    %add only offset errors and calculate stroke volume
    SV_MV_offsetonly(n) =  trapz(datan.MV.time, trueData.MV.mean + e.mv.systematic );%+ e.RR.MV
    SV_AV_offsetonly(n) =  trapz(datan.MV.time, trueData.AV.mean + e.av.systematic);% + e.RR.AV
    SV_AA_offsetonly(n) =  trapz(datan.MV.time, trueData.AA.mean + e.aa.systematic);% + e.RR.AA
    SV_PV_offsetonly(n) =  trapz(datan.MV.time, trueData.PV.mean + e.pv.systematic);% + e.RR.PV

    %add only systematic errors and calculate stroke volume
    SV_MV_systonly(n) =  trapz(datan.MV.time, trueData.MV.mean + e.mv.systematic + e.RR.MV);%
    SV_AV_systonly(n) =  trapz(datan.MV.time, trueData.AV.mean + e.av.systematic + e.RR.AV);% 
    SV_AA_systonly(n) =  trapz(datan.MV.time, trueData.AA.mean + e.aa.systematic + e.RR.AA);%
    SV_PV_systonly(n) =  trapz(datan.MV.time, trueData.PV.mean + e.pv.systematic + e.RR.PV);%
end

ctot_meanerr = mean(Ctoterr);
ctot_meanerr_percent = 100*mean(Ctoterr)/trueData.parameters.Ctot.mean;
ctot_sderr = std(Ctoterr);
Ctot_percentSD = 100*ctot_sderr/trueData.parameters.Ctot.mean;

rtot_meanerr = mean(Rtoterr);
rtot_meanerr_percent = 100*mean(Rtoterr)/trueData.parameters.Rtot.mean;
rtot_sderr = std(Rtoterr);
Rtot_percentSD = 100*rtot_sderr/trueData.parameters.Rtot.mean;

Emax_LVmeanerr = mean(Emax_LVerr);
Emax_LVmeanerr_percent = 100*mean(Emax_LVerr)/trueData.parameters.Emax_LV.mean;
Emax_LVsderr = std(Emax_LVerr);
Emax_percentSD = 100*Emax_LVsderr/trueData.parameters.Emax_LV.mean;

all=[ctot_meanerr,ctot_sderr,ctot_meanerr_percent,Ctot_percentSD;
    rtot_meanerr,rtot_sderr,rtot_meanerr_percent,Rtot_percentSD;
    Emax_LVmeanerr,Emax_LVsderr,Emax_LVmeanerr_percent,Emax_percentSD];
paramErrors = array2table(all,'Variablenames',{'mean','sd','mean (%)','sd (%)'},'rownames',{'Ctot','Rtot','Emax_LV'})

%% load a flow curve
% from true simulated data
load('dataSimulated.mat','dataSimulated')
trueData = dataSimulated.('simRandomSystematic');
origflowcurves = {trueData.MV.mean,trueData.AV.mean,trueData.AA.mean,trueData.PV.mean};
origflownames = {'MV','AV','AA','PV'};
plotnames = {'mitral valve','aortic valve','asc. aorta','pulm. veins'};
origTime = trueData.AV.time;

% "real" RR error
[~,e] = loadSampledMeasurementError(1,'simRandomSystematic');
realRRerrors =  {e.RR.MV,e.RR.AV,e.RR.AA,e.RR.PV};


%% Find mean of all sampled RR errors
allmeanRRerrors = cell(1,4);
realRRerrors = cell(1,100);
for f= 1:length(origflownames)
    for n = 1:100
        [datasampledn,e] = loadSampledMeasurementError(n,'simRandomSystematic');
        realRRerrors{n} =  {e.RR.MV,e.RR.AV,e.RR.AA,e.RR.PV};
        smoothingErrorReal = realRRerrors{n}{f};
        if n ==1
            allerrors = zeros(100,length(realRRerrors{n}{f}));
        end
        allerrors(n,:) = smoothingErrorReal';
    end
    %mean error
    meanRRerror = mean(allerrors);
    allmeanRRerrors{f} = meanRRerror';
end

realRRerrors = allmeanRRerrors;


%% extract rr data
for f = 1:length(origflowcurves)
    origFlow = origflowcurves{f};    
    % add temporal smoothing (RR variation) error
    rrdata.(origflownames{f}).mean = origFlow;
    rrdata.(origflownames{f}).time = origTime;
    rrdata.indtdiast = trueData.extra.indtdiast;
end

%% Create figure of +- sd for all errors
figure('Name','Fig4_Measurementerrors')
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 20 + 6;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(6+1,2,'TileSpacing','compact','Padding','compact')

letters = 'A':'F';
letternum = 1;

for f = 1:2
    [resolutionError,postPrError,smoothingErrorReal,offsetError,allErrors,origFlow,randerrors] = flowerror(f,origflowcurves,origflownames,realRRerrors,rrdata,origTime,trueData);
    nexttile
    hold on
    title(origflownames{f})
    errorbar(origTime,origFlow,allErrors,'.','color',[0.5 0.5 0.5],'LineWidth',1.2)
    errorbar(origTime,origFlow,randerrors,'k.','LineWidth',1.2)
    ylabel('Blood flow (mL/s)')
    set(gca,'FontSize',9)
    if strcmp('PV',origflownames{f})
        ylim([-40 180])
        yticks([-40 0 100 180])
    else
        ylim([-10 350])
        yticks([0 100 200 350])
    end
    xticks([origTime(1),round(origTime(end),1)])

    t=title(letters(letternum),'FontSize',12);
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';%TitleFontWeight
    pos = t.Position;
    pos(1)  = pos(1)-0.125;
    t.Position = pos;
    letternum = letternum+1;

end
legend({sprintf('Full uncertainty\n(for model estimation)'),'Random errors'}, ...
        'Position',[0.717 0.913404510596834 0.242528023814799 0.0599206362073385])

% EXAMPLE ERRORS
for f = 1:2
    [resolutionError,postPrError,smoothingErrorReal,offsetError,allErrors,origFlow,randerrors] = flowerror(f,origflowcurves,origflownames,realRRerrors,rrdata,origTime,trueData);
    % Plot all in one
    RRpos = smoothingErrorReal;
    RRpos(RRpos<0) = NaN;%0;
    RRneg = smoothingErrorReal;
    RRneg(RRneg>0) = NaN;%0;

    nexttile
    hold on
    neg = RRneg-offsetError-resolutionError-postPrError;
    neg(isnan(neg))=0;
    pos = RRpos+offsetError+resolutionError+postPrError;
    pos(isnan(pos))=0;
    plot(origTime, neg+pos,'-','color',[0.5 0.5 0],'LineWidth',2)

    fill([origTime;flip(origTime)],[offsetError+resolutionError+postPrError;flip(resolutionError+offsetError)],[0.8 0.8 0.8],'FaceAlpha',1,'Edgecolor','none')
    fill([origTime;flip(origTime)],[offsetError+resolutionError;flip(ones(size(origTime)).*offsetError)],[1 0.3 0.3],'FaceAlpha',1,'Edgecolor','none')
    fill([origTime;flip(origTime)],[ones(size(origTime)).*offsetError;zeros(size(origTime))],[0.5 0.5 1] ,'FaceAlpha',1,'Edgecolor','none')
    %neg:
    fill([origTime;flip(origTime)],[-offsetError-resolutionError-postPrError;flip(-resolutionError-offsetError)],[0.8 0.8 0.8],'FaceAlpha',1,'Edgecolor','none')
    fill([origTime;flip(origTime)],[-offsetError-resolutionError;flip(ones(size(origTime)).*-offsetError)],[1 0.3 0.3],'FaceAlpha',1,'Edgecolor','none')
    fill([origTime;flip(origTime)],[ones(size(origTime)).*-offsetError;zeros(size(origTime))],[0.5 0.5 1] ,'FaceAlpha',1,'Edgecolor','none')

    set(gca,'FontSize',9)
    ylim([-55 50])
    yticks([-55 0 25 50])
    xticks([origTime(1),round(origTime(end),1)])
    ylabel('Error (mL/s)')
    xlabel('Time (s)')
    yline(0);
end
legend({'RR variation','Observer variability','Spatial resolution','Offset'}, ...
            'Position',[0.735 0.74850784660891 0.223848073182635 0.0792328059484089])

%% RR variation example
%calc all heartbeats
doPlot = 0;
rrNumber=400;
randomDistribution = 1;
load('dataSimulated.mat','simulatedDataTrue');
fprintf('Simulating %d heartbeats...\n',rrNumber)
[sol,solTrue,allHeartbeats] = simulate_avatar_RR(simulatedDataTrue.allParameters,simulatedDataTrue.constants,simulatedDataTrue.options,simulatedDataTrue.ind,simulatedDataTrue.simtime,simulatedDataTrue.modelName,rrNumber,randomDistribution,doPlot);
Trange = allHeartbeats.t(:,end);
HRrange = 60./Trange;

[Trangesorted,sind] = sort(Trange);
HRrangesorted = HRrange(sind);
trueT = solTrue.t(end);
trueHR = 60/trueT;

col1 = [0.9 0.9 0.4];
col2 = [0.6 0.6 0.5];

hrcols = [linspace(col1(1),col2(1),length(Trange))',linspace(col1(2),col2(2),length(Trange))',linspace(col1(3),col2(3),length(Trange))'];

ind.MV = 4;
ind.AV = 6;
nexttile
hold on
for j = 1:40:length(Trange) %10 of the 400
    i = sind(j);
    tplots(j) = plot(allHeartbeats.t(i,:),allHeartbeats.x(i,:,ind.MV),'-','LineWidth',0.9,'color',hrcols(j,:));
end
trueplot = plot(solTrue.t,solTrue.x(:,ind.MV),'-','LineWidth',2,'color',[0 0 0]);
averageplot = plot(sol.t,sol.x(:,ind.MV),'--','LineWidth',2,'color',[0.5 0.5 0]);
xlabel('Time (s)')
ylabel('Blood flow (mL/s)')
xlim([0 max(Trange)])
ylim([-10 450])
yticks([0 100 200 350])

nexttile
hold on
for j = 1:40:length(Trange)
    i = sind(j);
    tplots(j) = plot(allHeartbeats.t(i,:),allHeartbeats.x(i,:,ind.AV),'-','LineWidth',0.9,'color',hrcols(j,:));
end
trueplot = plot(solTrue.t,solTrue.x(:,ind.AV),'-','LineWidth',2,'color',[0 0 0]);
averageplot = plot(sol.t,sol.x(:,ind.AV),'--','LineWidth',2,'color',[0.5 0.5 0]);
xlabel('Time (s)')
ylabel('Blood flow (mL/s)')
xlim([0 max(Trange)])
legend([tplots(1) tplots(end) trueplot averageplot],{sprintf('HR = %0.0f',HRrangesorted(1)),sprintf('HR = %0.0f',HRrangesorted(end)),sprintf('HR = %0.0f\n(True mean)',trueHR),'Averaged'},'Location','Northeast');
ylim([-10 350])
yticks([0 100 200 350])

%%
for f = 3:4
    [resolutionError,postPrError,smoothingErrorReal,offsetError,allErrors,origFlow,randerrors] = flowerror(f,origflowcurves,origflownames,realRRerrors,rrdata,origTime,trueData);
    nexttile
    hold on
    title(origflownames{f})
    errorbar(origTime,origFlow,allErrors,'.','color',[0.5 0.5 0.5],'LineWidth',1.2)
    errorbar(origTime,origFlow,randerrors,'k.','LineWidth',1.2)
    ylabel('Blood flow (mL/s)')
    set(gca,'FontSize',9)
    xticks([origTime(1),round(origTime(end),1)])
    if strcmp('PV',origflownames{f})
        ylim([-40 180])
        yticks([-40 0 100 180])
    else
        ylim([-10 350])
        yticks([0 100 200 350])
    end


    t=title(letters(letternum),'FontSize',12);
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    pos = t.Position;
    pos(1)  = pos(1)-0.125;
    t.Position = pos;
    letternum = letternum+1;
end

for f = 3:4
    [resolutionError,postPrError,smoothingErrorReal,offsetError,allErrors,origFlow,randerrors] = flowerror(f,origflowcurves,origflownames,realRRerrors,rrdata,origTime,trueData);
    RRpos = smoothingErrorReal;
    RRpos(RRpos<0) = NaN;%0;
    RRneg = smoothingErrorReal;
    RRneg(RRneg>0) = NaN;%0;

    % Plot all in one
    nexttile
    hold on

    neg = RRneg-offsetError-resolutionError-postPrError;
    neg(isnan(neg))=0;
    pos = RRpos+offsetError+resolutionError+postPrError;
    pos(isnan(pos))=0;
    plot(origTime, neg+pos,'-','color',[0.5 0.5 0],'LineWidth',2)

    fill([origTime;flip(origTime)],[offsetError+resolutionError+postPrError;flip(resolutionError+offsetError)],[0.8 0.8 0.8],'FaceAlpha',1,'Edgecolor','none')
    fill([origTime;flip(origTime)],[offsetError+resolutionError;flip(ones(size(origTime)).*offsetError)],[1 0.3 0.3],'FaceAlpha',1,'Edgecolor','none')
    fill([origTime;flip(origTime)],[ones(size(origTime)).*offsetError;zeros(size(origTime))],[0.5 0.5 1] ,'FaceAlpha',1,'Edgecolor','none')
    %neg:
    fill([origTime;flip(origTime)],[-offsetError-resolutionError-postPrError;flip(-resolutionError-offsetError)],[0.8 0.8 0.8],'FaceAlpha',1,'Edgecolor','none')
    fill([origTime;flip(origTime)],[-offsetError-resolutionError;flip(ones(size(origTime)).*-offsetError)],[1 0.3 0.3],'FaceAlpha',1,'Edgecolor','none')
    fill([origTime;flip(origTime)],[ones(size(origTime)).*-offsetError;zeros(size(origTime))],[0.5 0.5 1] ,'FaceAlpha',1,'Edgecolor','none')

    set(gca,'FontSize',9)
    ylim([-55 50])
    yticks([-55 0 25 50])
    xticks([origTime(1),round(origTime(end),1)])
    ylabel('Error (mL/s)')
    xlabel('Time (s)')
    xlabel('Time (s)')
    yline(0);
end

%% stroke volumes
nexttile([2 1])
hold on
ylabel('Stroke volume error (%)')
mverrAbs_syst = SV_MV_systonly-trapz(trueData.MV.time,trueData.MV.mean);
averrAbs_syst = SV_AV_systonly-trapz(trueData.AV.time,trueData.AV.mean);
aaerrAbs_syst = SV_AA_systonly-trapz(trueData.AA.time,trueData.AA.mean);
pverrAbs_syst = SV_PV_systonly-trapz(trueData.PV.time,trueData.PV.mean);
mverrAbs_systP = 100*mverrAbs_syst./trapz(trueData.MV.time,trueData.MV.mean);
averrAbs_systP = 100*averrAbs_syst./trapz(trueData.AV.time,trueData.AV.mean);
aaerrAbs_systP = 100*aaerrAbs_syst./trapz(trueData.AA.time,trueData.AA.mean);
pverrAbs_systP = 100*pverrAbs_syst./trapz(trueData.PV.time,trueData.PV.mean);
swarmchart(1.2.*ones(size(mverrAbs_systP)),mverrAbs_systP,2,[0.5 0.5 0],'*','XJitterWidth',0.10)
swarmchart(2.2.*ones(size(mverrAbs_systP)),averrAbs_systP,2,[0.5 0.5 0],'*','XJitterWidth',0.10)
swarmchart(3.2.*ones(size(mverrAbs_systP)),aaerrAbs_systP,2,[0.5 0.5 0],'*','XJitterWidth',0.10)
ss=swarmchart(4.2.*ones(size(mverrAbs_systP)),pverrAbs_systP,2,[0.5 0.5 0],'*','XJitterWidth',0.10);

mverr = 100*(SV_MV-trapz(trueData.MV.time,trueData.MV.mean))./trapz(trueData.MV.time,trueData.MV.mean);
averr = 100*(SV_AV-trapz(trueData.AV.time,trueData.AV.mean))./trapz(trueData.AV.time,trueData.AV.mean);
aaerr = 100*(SV_AA-trapz(trueData.AA.time,trueData.AA.mean))./trapz(trueData.AA.time,trueData.AA.mean);
pverr = 100*(SV_PV-trapz(trueData.PV.time,trueData.PV.mean))./trapz(trueData.PV.time,trueData.PV.mean);
s=swarmchart(1.*ones(size(mverr)),mverr,2,[0.5 0.5 1],'*','XJitterWidth',0.10);
swarmchart(2.*ones(size(mverr)),averr,2,[0.5 0.5 1],'*','XJitterWidth',0.10)
swarmchart(3.*ones(size(mverr)),aaerr,2,[0.5 0.5 1],'*','XJitterWidth',0.10)
swarmchart(4.*ones(size(mverr)),pverr,2,[0.5 0.5 1],'*','XJitterWidth',0.10)
yline(0,'--');

errorbar(1,mean(mverr),std(mverr),'k*','LineWidth',1,'markersize',5)
errorbar(2,mean(averr),std(averr),'k*','LineWidth',1,'markersize',5)
errorbar(3,mean(aaerr),std(aaerr),'k*','LineWidth',1,'markersize',5)
er=errorbar(4,mean(pverr),std(pverr),'k*','LineWidth',1,'markersize',5);

mverrAbs_rand = SV_MV_randonly-trapz(trueData.MV.time,trueData.MV.mean);
averrAbs_rand = SV_AV_randonly-trapz(trueData.AV.time,trueData.AV.mean);
aaerrAbs_rand = SV_AA_randonly-trapz(trueData.AA.time,trueData.AA.mean);
pverrAbs_rand = SV_PV_randonly-trapz(trueData.PV.time,trueData.PV.mean);
mverrAbs_randP = 100*mverrAbs_rand./trapz(trueData.MV.time,trueData.MV.mean);
averrAbs_randP = 100*averrAbs_rand./trapz(trueData.AV.time,trueData.AV.mean);
aaerrAbs_randP = 100*aaerrAbs_rand./trapz(trueData.AA.time,trueData.AA.mean);
pverrAbs_randP = 100*pverrAbs_rand./trapz(trueData.PV.time,trueData.PV.mean);
swarmchart(1.4.*ones(size(mverrAbs_randP)),mverrAbs_randP,2,[0.4 0 0.2],'*','XJitterWidth',0.10)
swarmchart(2.4.*ones(size(mverrAbs_randP)),averrAbs_randP,2,[0.4 0 0.2],'*','XJitterWidth',0.10)
swarmchart(3.4.*ones(size(mverrAbs_randP)),aaerrAbs_randP,2,[0.4 0 0.2],'*','XJitterWidth',0.10)
sr=swarmchart(4.4.*ones(size(mverrAbs_randP)),pverrAbs_randP,2,[0.4 0 0.2],'*','XJitterWidth',0.10);

xticks([1:4])
xticklabels({'MV','AV','AA','PV'})
yline(5,'--');
yline(-5,'--');
set(gca,'FontSize',9)
xlim([0.5 4.5])
t=title('E','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.6;
t.Position = pos;


nexttile([2 1])
hold on
ylabel('Stroke volume difference (%)')
diffMVAV_systonly = 100*(SV_MV_systonly-SV_AV_systonly)./mean([SV_MV_systonly,SV_AV_systonly],2);
diffMVAA_systonly = 100*(SV_MV_systonly-SV_AA_systonly)./mean([SV_MV_systonly,SV_AA_systonly],2);
diffMVPV_systonly = 100*(SV_MV_systonly-SV_PV_systonly)./mean([SV_MV_systonly,SV_PV_systonly],2);
diffAVAA_systonly = 100*(SV_AV_systonly-SV_AA_systonly)./mean([SV_AV_systonly,SV_AA_systonly],2);
diffAVPV_systonly = 100*(SV_AV_systonly-SV_PV_systonly)./mean([SV_AV_systonly,SV_PV_systonly],2);
diffAAPV_systonly = 100*(SV_AA_systonly-SV_PV_systonly)./mean([SV_AA_systonly,SV_PV_systonly],2);
swarmchart(1.2.*ones(size(diffMVAV_systonly)),diffMVAV_systonly,2,[0.5 0.5 0],'*','XJitterWidth',0.10)
swarmchart(2.2.*ones(size(diffMVAV_systonly)),diffMVAA_systonly,2,[0.5 0.5 0],'*','XJitterWidth',0.10)
swarmchart(3.2.*ones(size(diffMVAV_systonly)),diffMVPV_systonly,2,[0.5 0.5 0],'*','XJitterWidth',0.10)
ss=swarmchart(4.2.*ones(size(diffMVAV_systonly)),diffAVAA_systonly,2,[0.5 0.5 0],'*','XJitterWidth',0.10);
swarmchart(5.2.*ones(size(diffMVAV_systonly)),diffAVPV_systonly,2,[0.5 0.5 0],'*','XJitterWidth',0.10)
swarmchart(6.2.*ones(size(diffMVAV_systonly)),diffAAPV_systonly,2,[0.5 0.5 0],'*','XJitterWidth',0.10)

diffMVAV = 100*(SV_MV-SV_AV)./mean([SV_MV,SV_AV],2);
diffMVAA = 100*(SV_MV-SV_AA)./mean([SV_MV,SV_AA],2);
diffMVPV = 100*(SV_MV-SV_PV)./mean([SV_MV,SV_PV],2);
diffAVAA = 100*(SV_AV-SV_AA)./mean([SV_AV,SV_AA],2);
diffAVPV = 100*(SV_AV-SV_PV)./mean([SV_AV,SV_PV],2);
diffAAPV = 100*(SV_AA-SV_PV)./mean([SV_AA,SV_PV],2);

swarmchart(1.*ones(size(SV_MV)),100*(SV_MV-SV_AV)./SV_MV,2,[0.5 0.5 1],'*','XJitterWidth',0.10)
swarmchart(2.*ones(size(SV_MV)),100*(SV_MV-SV_AA)./SV_MV,2,[0.5 0.5 1],'*','XJitterWidth',0.10)
swarmchart(3.*ones(size(SV_MV)),100*(SV_MV-SV_PV)./SV_MV,2,[0.5 0.5 1],'*','XJitterWidth',0.10)
swarmchart(4.*ones(size(SV_MV)), 100*(SV_AV-SV_AA)./SV_AV,2,[0.5 0.5 1],'*','XJitterWidth',0.10)
swarmchart(5.*ones(size(SV_MV)), 100*(SV_AV-SV_PV)./SV_AV,2,[0.5 0.5 1],'*','XJitterWidth',0.10)
swarmchart(6.*ones(size(SV_MV)), 100*(SV_AA-SV_PV)./SV_AA,2,[0.5 0.5 1],'*','XJitterWidth',0.10)

errorbar(1,mean(diffMVAV),std(diffMVAV),'k*','LineWidth',1,'markersize',5)
errorbar(2,mean(diffMVAA),std(diffMVAA),'k*','LineWidth',1,'markersize',5)
errorbar(3,mean(diffMVPV),std(diffMVPV),'k*','LineWidth',1,'markersize',5)
errorbar(4,mean(diffAVAA),std(diffAVAA),'k*','LineWidth',1,'markersize',5)
errorbar(5,mean(diffAVPV),std(diffAVPV),'k*','LineWidth',1,'markersize',5)
errorbar(6,mean(diffAAPV),std(diffAAPV),'k*','LineWidth',1,'markersize',5)

diffMVAV_randonly = 100*(SV_MV_randonly-SV_AV_randonly)./mean([SV_MV_randonly,SV_AV_randonly],2);
diffMVAA_randonly = 100*(SV_MV_randonly-SV_AA_randonly)./mean([SV_MV_randonly,SV_AA_randonly],2);
diffMVPV_randonly = 100*(SV_MV_randonly-SV_PV_randonly)./mean([SV_MV_randonly,SV_PV_randonly],2);
diffAVAA_randonly = 100*(SV_AV_randonly-SV_AA_randonly)./mean([SV_AV_randonly,SV_AA_randonly],2);
diffAVPV_randonly = 100*(SV_AV_randonly-SV_PV_randonly)./mean([SV_AV_randonly,SV_PV_randonly],2);
diffAAPV_randonly = 100*(SV_AA_randonly-SV_PV_randonly)./mean([SV_AA_randonly,SV_PV_randonly],2);
swarmchart(1.4.*ones(size(diffMVAV_randonly)),diffMVAV_randonly,2,[0.4 0 0.2],'*','XJitterWidth',0.10)
swarmchart(2.4.*ones(size(diffMVAV_randonly)),diffMVAA_randonly,2,[0.4 0 0.2],'*','XJitterWidth',0.10)
swarmchart(3.4.*ones(size(diffMVAV_randonly)),diffMVPV_randonly,2,[0.4 0 0.2],'*','XJitterWidth',0.10)
sr=swarmchart(4.4.*ones(size(diffMVAV_randonly)),diffAVAA_randonly,2,[0.4 0 0.2],'*','XJitterWidth',0.10);
swarmchart(5.4.*ones(size(diffMVAV_randonly)),diffAVPV_randonly,2,[0.4 0 0.2],'*','XJitterWidth',0.10)
swarmchart(6.4.*ones(size(diffMVAV_randonly)),diffAAPV_randonly,2,[0.4 0 0.2],'*','XJitterWidth',0.10)

yline(5,'--');
yline(-5,'--');
yline(0,'--');

xlim([0.5 6.5])
xticks([1:6])
xticklabels({'MV-AV','MV-AA','MV-PV','AV-AA','AV-PV','AA-PV'})
legend([s(1);er;sr(1);ss(1)],{'All errors','Mean+-sd of all','Random error',...
    'Systematic error'},'location','northoutside','Numcolumns',2,...
    'Position',[0.0872432573356074 0.284542914516151 0.4 0.04])
set(gca,'FontSize',9)
t=title('F','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.8;
t.Position = pos;




%% SV tables
mverrAbs = SV_MV-trapz(trueData.MV.time,trueData.MV.mean);
averrAbs = SV_AV-trapz(trueData.AV.time,trueData.AV.mean);
aaerrAbs = SV_AA-trapz(trueData.AA.time,trueData.AA.mean);
pverrAbs = SV_PV-trapz(trueData.PV.time,trueData.PV.mean);

meanandsds =[mean(mverr),std(mverr),mean(mverrAbs),std(mverrAbs);...
    mean(averr),std(averr),mean(averrAbs),std(averrAbs);...
    mean(aaerr),std(aaerr),mean(aaerrAbs),std(aaerrAbs);...
    mean(pverr),std(pverr),mean(pverrAbs),std(pverrAbs)];
SVerrorsTable = array2table(meanandsds,'Rownames',{'MV','AV','AA','PV'},'VariableNames',{'Mean error (%)','SD error (%)', 'Mean error (abs)','SD error (abs)'})

internaldiffs = [std(SV_MV-SV_AV),mean(diffMVAV),std(diffMVAV),std(diffMVAV_randonly),std(diffMVAV_systonly);...
    std(SV_MV-SV_AA), mean(diffMVAA),std(diffMVAA),std(diffMVAA_randonly),std(diffMVAA_systonly);...
    std(SV_MV-SV_PV),mean(diffMVPV),std(diffMVPV),std(diffMVPV_randonly),std(diffMVPV_systonly);...
    std(SV_AV-SV_AA), mean(diffAVAA),std(diffAVAA),std(diffAVAA_randonly),std(diffAVAA_systonly);...
    std(SV_AV-SV_PV),mean(diffAVPV),std(diffAVPV),std(diffAVPV_randonly),std(diffAVPV_systonly);...
    std(SV_AA-SV_PV),mean(diffAAPV),std(diffAAPV),std(diffAAPV_randonly),std(diffAAPV_systonly)];
SVdiffsTable = array2table(internaldiffs,'Rownames',{'MV-AV','MV-AA','MV-PV','AV-AA','AV-PV','AA-PV'},'VariableNames',{'SD Abs diff (mL)','Mean diff (%)', 'sd diff (%)','sd diff random only (%)','sd diff syst only (%)'})



meanandsdsSystRand =[mean(mverrAbs_randP),std(mverrAbs_randP),mean(mverrAbs_systP),std(mverrAbs_systP);...
    mean(averrAbs_randP),std(averrAbs_randP),mean(averrAbs_systP),std(averrAbs_systP);...
    mean(aaerrAbs_randP),std(aaerrAbs_randP),mean(aaerrAbs_systP),std(aaerrAbs_systP);...
    mean(pverrAbs_randP),std(pverrAbs_randP),mean(pverrAbs_systP),std(pverrAbs_systP)];
SVerrorsRandSystTable = array2table(meanandsdsSystRand,'Rownames',{'MV','AV','AA','PV'},...
    'VariableNames',{'random Mean error (%)','random SD error (%)', 'systematic Mean error (%)','systematic sd error (%)'})


randerr = split(sprintf('%0.2f +- %0.2fx',meanandsdsSystRand(:,1:2)'),'x'); %in percent
systerr = split(sprintf('%0.2f +- %0.2fx',meanandsdsSystRand(:,3:4)'),'x'); %in percent
toterrAbs = split(sprintf('%0.2f +- %0.2fx',meanandsds(:,3:4)'),'x'); %abs
toterrPerc = split(sprintf('%0.2f +- %0.2fx',meanandsds(:,1:2)'),'x'); %percent
internaldiffs = [mean(SV_MV-SV_AV),std(SV_MV-SV_AV),mean(diffMVAV),std(diffMVAV),mean(diffMVAV_randonly),std(diffMVAV_randonly),mean(diffMVAV_systonly),std(diffMVAV_systonly);...
    mean(SV_MV-SV_AA),std(SV_MV-SV_AA), mean(diffMVAA),std(diffMVAA),mean(diffMVAA_randonly),std(diffMVAA_randonly),mean(diffMVAA_systonly),std(diffMVAA_systonly);...
    mean(SV_MV-SV_PV),std(SV_MV-SV_PV),mean(diffMVPV),std(diffMVPV),mean(diffMVPV_randonly),std(diffMVPV_randonly),mean(diffMVPV_systonly),std(diffMVPV_systonly);...
    mean(SV_AV-SV_AA),std(SV_AV-SV_AA), mean(diffAVAA),std(diffAVAA),mean(diffAVAA_randonly),std(diffAVAA_randonly),mean(diffAVAA_systonly),std(diffAVAA_systonly);...
    mean(SV_AV-SV_PV),std(SV_AV-SV_PV),mean(diffAVPV),std(diffAVPV),mean(diffAVPV_randonly),std(diffAVPV_randonly),mean(diffAVPV_systonly),std(diffAVPV_systonly);...
    mean(SV_AA-SV_PV),std(SV_AA-SV_PV),mean(diffAAPV),std(diffAAPV),mean(diffAAPV_randonly),std(diffAAPV_randonly),mean(diffAAPV_systonly),std(diffAAPV_systonly)];
diffsTotAbs = split(sprintf('%0.2f +- %0.2fx',internaldiffs(:,1:2)'),'x'); %abs diff total
diffsTotPerc = split(sprintf('%0.2f +- %0.2fx',internaldiffs(:,3:4)'),'x'); %perc diff total
diffsRand = split(sprintf('%0.2f +- %0.2fx',internaldiffs(:,5:6)'),'x'); %rand diff perc
diffsSyst = split(sprintf('%0.2f +- %0.2fx',internaldiffs(:,7:8)'),'x'); %syst diff perc

SVerrorsTableTEXT = table([toterrAbs(1:end-1);diffsTotAbs(1:end-1)],...
    [toterrPerc(1:end-1);diffsTotPerc(1:end-1)],...
    [randerr(1:end-1);diffsRand(1:end-1)], ...
    [systerr(1:end-1);diffsSyst(1:end-1)], ...
    'VariableNames',{'Total error (mL)','Total error (%)','Random error only (%)','Systematic error only (%)'}, ...
    'Rownames',{'MV','AV','AA','PV','MV-AV','MV-AA','MV-PV','AV-AA','AV-PV','AA-PV'})


%% save stuff
if saveStuff
    if nargin <3
        plotFolderName ='C:\Users\kajtu\Modeller\Uncertainty-estimation\Plots\Dataplots\SampledDataErrors';
    end
    saveAllFigures(plotFolderName,saveStuff,1)
    writetable(paramErrors,fullfile(plotFolderName,'paramErrors.xlsx'),"WriteRowNames",1)
    writetable(SVerrorsTableTEXT,fullfile(plotFolderName,'SVerrorsTableTEXT.xlsx'),"WriteRowNames",1)
end

end



function [resolutionError,postPrError,smoothingErrorReal,offsetError,allErrors,origFlow,randerrors] = flowerror(f,origflowcurves,origflownames,realRRerrors,rrdata,origTime,trueData)
    origFlow = origflowcurves{f};
    % add 0.9 cm/s * vessel area offset
    if strcmp('MV',origflownames{f}) || strcmp('AA',origflownames{f})
        offsetError =  0.9*7;
    elseif strcmp('AV',origflownames{f})
        offsetError =  0.9*5;
    elseif strcmp('PV',origflownames{f})
        offsetError = 0.9*5.3;
    end

    % add spatial resolution error
    resolutionErrorPercent = 0.02;% 2%
    if strcmp(origflownames{f} , 'PV')
        resolutionErrorPercent = resolutionErrorPercent*2;% 4%
    end
    resolutionError = abs(origFlow).*resolutionErrorPercent;

    % add post processing error (Observer variability in post processing)
    postProcessingErrorPercent = 0.05; % 5%
    if strcmp(origflownames{f} , 'PV')
        postProcessingErrorPercent = 0.07;% +2% for the more complicated method and summation of 3-4 veins
    end
    postPrError = abs(origFlow).*postProcessingErrorPercent;

    % add temporal smoothing (RR variation) error
    smoothingErrorReal = realRRerrors{f};

    [smoothingError] = calcRRsmoothingerror(rrdata);
    smoothingError = smoothingError.(origflownames{f}).mean;

    % random errors  
    randompercent = sqrt(resolutionErrorPercent^2 + postProcessingErrorPercent^2);
    randerrors = abs(origFlow).*randompercent;
    randerrors(randerrors < mean(randerrors)) = mean(randerrors);

    % all errors in one (to use as errorbars)
    allErrors = randerrors + offsetError + abs(smoothingError);
end
function plotSampledData(experimentName,saveFigures,plotFolderName)
% sample bootstrapped data and plot it

if nargin <1
    experimentName = 'simRandomSystematic';%simRandomSystematic or simRandomOnly
end
if nargin < 2
        saveFigures =1;
end

load('dataSimulated.mat','dataSimulated')
trueData = dataSimulated.(experimentName);

estimationData = cell(1,10);
e=1;
cd ../Optimization
for i = 1:10:100 % 10 selected samples
    try
    [estimationData{e},~] = loadSampledMeasurementError(i,experimentName);
    catch
        disp(['no more sampled data found (i = ' num2str(i)])
        estimationData(e:end) = [];
        break
    end
    e = e+1;
end
cd ../Plotting

num2plot = 10;
step = round(length(estimationData)/num2plot);

experimentcolors = tab10(length(estimationData));
colors = experimentcolors;
greenCol = [0 0 0];

%% Create figure
figure('Name',sprintf('True simulated data + %d sampled (%s)',length(estimationData),experimentName))
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 17;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(3,6,'TileSpacing','loose','Padding','compact')

nexttile([1,3])
hold on
for i = 1:length(estimationData)
        plot(estimationData{i}.MV.time,estimationData{i}.MV.mean,'.-','color',colors(i,:),'Markersize',12,'LineWidth',1)
end
plot(trueData.MV.time,trueData.MV.mean,'o','color',greenCol,'linewidth',2,'Markersize',5)
ylabel('Blood flow in MV (mL/s)')
xlabel('Time (s)')
tend = trueData.MV.time(end);
xlim([0 tend])
set(gca,'FontSize',10)
t=title('A','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.16;
t.Position = pos;

nexttile([1,3])
hold on
for i = 1:length(estimationData)
        est(i)=plot(estimationData{i}.AV.time,estimationData{i}.AV.mean,'.-','color',colors(i,:),'Markersize',12,'LineWidth',1);
end
trueplot = plot(trueData.AV.time,trueData.AV.mean,'o','color',greenCol,'linewidth',2,'Markersize',5);
ylabel('Blood flow in AV (mL/s)')
xlabel('Time (s)')
n = split(sprintf('Sampled data %dx',1:length(estimationData)),'x');
legend([est,trueplot],[n(1:end-1);{'True data'}], ...
    'Position',[0.75395269610321 0.754276960784314 0.22 0.23])

xlim([0 tend])
set(gca,'FontSize',10)
t=title('B','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.16;
t.Position = pos;

nexttile([1,3])
hold on
for i = 1:length(estimationData)
        plot(estimationData{i}.AA.time,estimationData{i}.AA.mean,'.-','color',colors(i,:),'Markersize',12,'LineWidth',1)
end
plot(trueData.AA.time,trueData.AA.mean,'o','color',greenCol,'linewidth',2,'Markersize',5)
ylabel('Blood flow in AA (mL/s)')
xlabel('Time (s)')
xlim([0 tend])
set(gca,'FontSize',10)
t=title('C','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.16;
t.Position = pos;

nexttile([1,3])
hold on
for i = 1:length(estimationData)
        plot(estimationData{i}.PV.time,estimationData{i}.PV.mean,'.-','color',colors(i,:),'Markersize',12,'LineWidth',1)
end
plot(trueData.PV.time,trueData.PV.mean,'o','color',greenCol,'linewidth',2,'Markersize',5)
ylabel('Blood flow in PV (mL/s)')
xlabel('Time (s)')
xlim([0 tend])
set(gca,'FontSize',10)
t=title('D','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.16;
t.Position = pos;


nexttile
hold on
for i = 1:step:length(estimationData)
    plot(1,estimationData{i}.SBP.mean,'.','color',colors(i,:),'linewidth',0.9,'Markersize',14)
    plot(2,estimationData{i}.DBP.mean,'.','color',colors(i,:),'linewidth',0.9,'Markersize',14)
end
plot(1,trueData.SBP.mean,'o','color',greenCol,'linewidth',2,'Markersize',5)
plot(2,trueData.DBP.mean,'o','color',greenCol,'linewidth',2,'Markersize',5)
ylabel('Blood pressure (mmHg)')
xticks([1,2])
xticklabels({'SBP','DBP'})
xlim([0 3])
set(gca,'FontSize',10)
t=title('E','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.31;
t.Position = pos;

pnames = fieldnames(trueData.parameters);
letters = 'F':'Z';
for p = 1:length(pnames)
    nexttile
    hold on
    for i = 1:length(estimationData)
        plot(p,estimationData{i}.parameters.(pnames{p}).mean,'.','color',colors(i,:),'linewidth',1,'Markersize',14)
    end
    plot(p,trueData.parameters.(pnames{p}).mean,'o','color',greenCol,'linewidth',2,'Markersize',5)
    ylabel(pnames{p})
    xticks([])
    set(gca,'FontSize',10)
    t=title(letters(p),'FontSize',12);
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    pos = t.Position;
    pos(1)  = pos(1)-1.5;
    t.Position = pos;
end

%% save plot
if saveFigures
    if nargin >2
        plotFolder = plotFolderName;
    else
        plotFolderbase = 'C:\Users\kajtu\Modeller\Uncertainty-estimation\Plots\Dataplots';
        plotFolder = fullfile(plotFolderbase,[experimentName '_' datestr(now,'yymmdd-HHMM')]);
        mkdir(plotFolder)
    end
    saveAllFigures(plotFolder,1,1)
end

end

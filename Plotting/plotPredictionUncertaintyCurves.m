function plotPredictionUncertaintyCurves(simulations,ynames,figureName,colors,trueData)

plotTrueData = nargin > 4;
fz1 = 8.5;
fz2 = 10;
green = [0 0.55 0.35];
grey = green;
yellow = [0.6 0.4 0];

r = 3;c=2;
figure('Name',figureName,'Visible','on')
set(gcf,'Color','white')
set(gcf, 'InvertHardCopy', 'off'); % setting 'grid color reset' off
xdim_CM = 17;
ydim_CM = 12.75;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiles=tiledlayout(r,c,'TileSpacing','loose','Padding','compact');

nexttile
yInd = strcmp('mvCorr',ynames);
hold on
for s = 1:length(simulations)
    for p = 1:length(simulations{s}.time)
        if length(simulations) == 1
            colind = p;
        else
            colind = s;
        end
        tend = simulations{s}.time{p}(end);
        t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
        yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
        fill(t,yvals,colors{colind},'FaceAlpha',0.4,'EdgeColor','none');
        plot(100.*simulations{s}.time{p}./tend,simulations{s}.Observables.best{p}(:,yInd),'-','color',colors{colind},'LineWidth',2.2)
    end
end
if plotTrueData
    errorbar(100.*trueData.MV.time./trueData.MV.time(end),trueData.MV.mean,trueData.MV.sem,'o','color',grey,'LineWidth',1.5,'MarkerSize',5)
end
xlabel('% of cardiac cycle','FontSize',fz1)
ylabel('Flow in mitral valve (ml/s)','FontSize',fz1)
set(gca,'FontSize',fz1)
t=title('A','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.2;
t.Position = pos;

nexttile
hold on
yInd = strcmp('avCorr',ynames);
for s = 1:length(simulations)
    for p = 1:length(simulations{s}.time)
        tend = simulations{s}.time{p}(end);
        t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
        yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
        if length(simulations) == 1
            colind = p;
            sims(p)=fill(t,yvals,colors{colind},'FaceAlpha',0.4,'EdgeColor','none');
        else
            colind = s;
            sims(s)=fill(t,yvals,colors{colind},'FaceAlpha',0.4,'EdgeColor','none');
        end
        simsbest = plot(100.*simulations{s}.time{p}./tend,simulations{s}.Observables.best{p}(:,yInd),'-','color',colors{colind},'LineWidth',2.2);
    end
end
if plotTrueData
    truep = errorbar(100.*trueData.AV.time./trueData.AV.time(end),trueData.AV.mean,trueData.AV.sem,'o','color',grey,'LineWidth',1.5,'MarkerSize',5);
end
xlabel('% of cardiac cycle','FontSize',fz1)
ylabel('Flow in aortic valve (ml/s)','FontSize',fz1)
set(gca,'FontSize',fz1)
t=title('B','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.2;
t.Position = pos;

t1=nexttile;
colororder(t1,{'k','k'})%y axis colors
yInd = strcmp('Ela',ynames);
hold on
yyaxis right
ymax = 0;
for s = 1:length(simulations)
    for p = 1:length(simulations{s}.time)
        if length(simulations) == 1;colind = p;
        else; colind = s;end
        tend = simulations{s}.time{p}(end);
        t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
        yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
        fill(t,yvals,colors{colind},'FaceAlpha',0.4,'EdgeColor','none')
        ymax = max([ymax,max(simulations{s}.Observables.max{p}(:,yInd))]);
    end
end
if plotTrueData
    truep2 = plot(100.*trueData.allSimulations.t./trueData.allSimulations.t(end),trueData.allSimulations.y(:,yInd),'--','color',grey,'LineWidth',2,'MarkerSize',6);
end
ylim([-ymax ymax])

yyaxis left
yInd = strcmp('Elv',ynames);
ymax=0;
for s = 1:length(simulations)
    for p = 1:length(simulations{s}.time)
        if length(simulations) == 1;colind = p;
        else; colind = s;end
        tend = simulations{s}.time{p}(end);
        t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
        yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
        fill(t,yvals,colors{colind},'FaceAlpha',0.4,'EdgeColor','none')
        ymax = max([ymax,max(simulations{s}.Observables.max{p}(:,yInd))]);
    end
end
if plotTrueData
    truep2 = plot(100.*trueData.allSimulations.t./trueData.allSimulations.t(end),trueData.allSimulations.y(:,yInd),'--','color',grey,'LineWidth',2,'MarkerSize',6);
end
ylim([0 ymax*2])
xlabel('% of cardiac cycle','FontSize',fz1)
yyaxis left
ylabel(sprintf('Time-varying elastance\nin LV (mmHg/ml)'),'FontSize',fz1)
yyaxis right
ylabel(sprintf('Time-varying elastance\nin LA (mmHg/ml)'),'FontSize',fz1)
set(gca,'FontSize',fz1)
xlim([0 100])
t=title('C','FontSize',12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';
pos = t.Position;
pos(1)  = pos(1)-0.2;
t.Position = pos;

yparamsToPlot = {'P_Aortic','pLA','pLV'};
ylabels = {'Aortic pressure (mmHg)','Pressure in LA (mmHg)','Pressure in LV (mmHg)'};
letters = 'D':'Z';
for y = 1:length(yparamsToPlot)
    nexttile
    yInd = strcmp(yparamsToPlot{y},ynames);
    hold on
    for s = 1:length(simulations)
        for p = 1:length(simulations{s}.time)
            if length(simulations) == 1;colind = p;
            else; colind = s;end
            tend = simulations{s}.time{p}(end);
            t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
            yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
            fill(t,yvals,colors{colind},'FaceAlpha',0.4,'EdgeColor','none')%0.6
        end
    end
    if plotTrueData
        truep = plot(100.*trueData.allSimulations.t./trueData.allSimulations.t(end),trueData.allSimulations.y(:,yInd),'--','color',grey,'LineWidth',2,'MarkerSize',6);
    end
    xlabel('% of cardiac cycle','FontSize',fz1)
    ylabel(ylabels{y},'FontSize',fz1)
    set(gca,'FontSize',fz1)
    xlim([0 100])
    t=title(letters(y),'FontSize',12);
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';%TitleFontWeight
    pos = t.Position;
    pos(1)  = pos(1)-0.2;
    t.Position = pos;
end

if plotTrueData
    l=legend([sims(1),simsbest(1), truep,truep2],{'Model uncertainty','Best fit','True values','True values'},'Location','NorthEast','FontSize',8);
end
end
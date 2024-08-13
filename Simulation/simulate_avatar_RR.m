function [sol,solTrue,allHeartbeats] = simulate_avatar_RR(theta,constants,options,ind,simtimeOrig,modelName,rrNumber,randomDistribution,doPlot,Trange)
% Function to simulate several heartbeats and averaging the resulting flow curves.
% rrNumber: number of heartbeats each simulation should average over
% simtimeOrig: timepoints in seconds for the original heart rate
% randomDistribution: if set to 1, the range of cardiac cycle lengths T is
% randomly sampled from a normal distribution. If not, a fixed normal
% distribution will be used.
% Trange: provided if a specific range of Ts is wanted, otherwise do not
% use this argument

dispErrors=1;
T = constants(ind.T);

% simulate the "true" heartrate
warning('off')
[solTrue] = simulate_avatar_ss(theta,constants,options,ind,simtimeOrig,modelName);
warning('on')


if nargin <10 %if no Trange provided
    %% gaussian distribution of Ts
    meanT = T;
    sdT = 0.07*T;

    if randomDistribution
        Trange = normrnd(meanT,sdT,[1 rrNumber]); % sample randomly from normal distibution
    else
        rrstep = 1/(rrNumber+1);
        ps = rrstep:rrstep:1-rrstep; %sample with a fixed distance from the normal distribution for consistency
        Trange=norminv(ps,meanT,sdT);
    end
end

%% simulate several heartbeats
allHeartbeats.x = zeros([length(Trange),size(solTrue.x)]);
allHeartbeats.y = zeros([length(Trange),size(solTrue.y)]);
allHeartbeats.normalizedtimes = zeros(length(Trange),length(simtimeOrig));
allHeartbeats.t = zeros(length(Trange),length(simtimeOrig));
percentSimtime = simtimeOrig./T;
failedSimulations = zeros(size(Trange));
k_diast_LV = zeros(size(Trange));
k_diast_LA = zeros(size(Trange));
Emax_LV = zeros(size(Trange));
k_syst_LA = zeros(size(Trange));
for i = 1:length(Trange)
    % Set T and simulation time
    simtime = percentSimtime .* Trange(i);
    constants(ind.T) = Trange(i);
    paramValues = theta;

    %Set kdiast LV and LA relative T = T at rest
    k=-0.55;
    k_diast_LV(i) =  Trange(i)*k + (theta(ind.k_diast_LV)-(T*k));
    paramValues(ind.k_diast_LV) = k_diast_LV(i);

    k=0; %decreases slightly with hr = increases slightly with T increase
    k_syst_LA(i) =  Trange(i)*k + (theta(ind.k_syst_LA)-(T*k));%  
    paramValues(ind.k_syst_LA) =  k_syst_LA(i);

    k_diast_LA(i) = theta(ind.k_diast_LA);
    paramValues(ind.k_diast_LA) =  k_diast_LA(i);

    %Emax LV larger in higher HR
    k = -6;
    Emax_LV(i) =  Trange(i)*k + (theta(ind.Emax_LV)-(T*k));
    paramValues(ind.Emax_LV) =  Emax_LV(i);

    %onset_LA
    k = 0.1;
    onset_LA(i) =  Trange(i)*k + (theta(ind.onset_LA)-(T*k));
    paramValues(ind.onset_LA) =  onset_LA(i);

    % simulate
    warning('off')
    sim = simulate_avatar_ss(paramValues,constants,options,ind,simtime,modelName);
    warning('on')
    nansum = sum(isnan(sim.x(:)));
    if (sim.status<0) || (nansum > 0) %if simulation is NaN (ie the simulation failed)
        %disp(['simulate_avatar_RR: Simulation failed. Nan sum ' num2str(nansum) ', Status ' num2str(sim.status)])
        failedSimulations(i) = 1;
    else
        allHeartbeats.x(i,:,:) = sim.x;
        allHeartbeats.y(i,:,:) = sim.y;
        allHeartbeats.t(i,:) = sim.t;
    end

    allHeartbeats.normalizedtimes(i,:) = sim.t./Trange(i);

end

%remove any NaN simulations
if sum(failedSimulations)>0
    if dispErrors
    fprintf('%d of %d simulations failed\n',sum(failedSimulations),length(Trange))
    disp(find(failedSimulations))
    end
    if sum(failedSimulations)/length(Trange) < 0.2 %if only a few failed, just remove them from the averaging
        failedSimulations = logical(failedSimulations);
        allHeartbeats.x(failedSimulations,:,:) = [];
        allHeartbeats.y(failedSimulations,:,:) = [];
        allHeartbeats.t(failedSimulations,:) = [];
        allHeartbeats.normalizedtimes(failedSimulations,:) = [];
        Trange(failedSimulations) = [];
    end
end

% average all heartbeats
sol = sim;
sol.t = simtimeOrig';
meanx = mean(allHeartbeats.x);
sol.x = squeeze(meanx);
meany = mean(allHeartbeats.y);
sol.y = squeeze(meany);

sdx = std(allHeartbeats.x);
sol.sdx = squeeze(sdx);
sdy = std(allHeartbeats.y);
sol.sdy = squeeze(sdy);

%% Plot 
if doPlot
        HRrange = 60./Trange;
        colors = magma(length(HRrange));
    if nargin <10
        x = 0.8*T:0.0001:1.2*T;
        y = normpdf(x,meanT,sdT);
        Trangerandom = normrnd(meanT,sdT,[1 rrNumber]); %if sample them randomly


        figure('Name','Trange_distribution_rrsim')
        set(gcf,'Color','white')
        xdim_CM = 10;
        ydim_CM = 15;
        set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
        set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
        tiledlayout('flow')

        nexttile
        plot(x,y,'-')
        hold on
        histogram(Trange)
        plot(Trangerandom,ones(size(Trangerandom)),'k*')
        plot(Trange,ones(size(Trange)),'go')
        xline(meanT,'--')
        xline(meanT+sdT,'--')
        xline(meanT-sdT,'--')
        xline(Trange,'-','color',[0.7 0.7 0.7])
        xlabel('T (s)')
        xlim([min(x),max(x)])

        nexttile
        hold on
        histogram(HRrange)
        plot(HRrange,ones(size(HRrange)),'go')
        xline(mean(HRrange),'--')
        xline(mean(HRrange)+std(HRrange),'--')
        xline(mean(HRrange)-std(HRrange),'--')
        xline(HRrange,'-','color',[0.7 0.7 0.7])
        xlabel('HR (beats/minute)')
    end

    figure('Name','True mean T vs smoothed_rrsim')
    toplot = {'MV','AV','AC','PV'};
    for i = 1:length(toplot)
        nexttile
        hold on
        plotind = ind.(toplot{i});
        plot(solTrue.t,solTrue.x(:,plotind),'k-','LineWidth',1.3)
        plot(sol.t,sol.x(:,plotind),'--','Color',[1 0.5 0.5],'LineWidth',1.3)
        xlabel('Time (s)')
        ylabel(['Blood flow' toplot{i} ' (mL/s)'])
    end

    figure('Name','Parameter change')
    nexttile
    hold on
    ylabel('k diast (*) &  k syst LA (^)')
    xlabel('Heart rate (beats/min)')
    for i = 1:length(Trange)
        plot(60/Trange(i),k_diast_LA(i),'*','color',colors(i,:))
        plot(60/Trange(i),k_syst_LA(i),'^','color',colors(i,:))
    end

    nexttile
    hold on
    ylabel('k diast LV')
    xlabel('Heart rate (beats/min)')
    for i = 1:length(Trange)
        plot(60/Trange(i),k_diast_LV(i),'*','color',colors(i,:))
    end

    nexttile
    hold on
    ylabel('Emax LV')
    xlabel('Heart rate (beats/min)')
    for i = 1:length(Trange)
        plot(60/Trange(i),Emax_LV(i),'*','color',colors(i,:))
    end
    
    % in percent
    nexttile
    xlabel('% of cardiac cycle')
    ylabel('MV')
    hold on
    for i = 1:length(Trange)
        plot(100.*simtimeOrig./T,allHeartbeats.x(i,:,ind.MV),'color',colors(i,:))
    end
    plot(100*(solTrue.t)./T,solTrue.x(:,ind.MV),'k-','LineWidth',1.3)
    plot(100*sol.t./sol.t(end),sol.x(:,ind.MV),'--','Color',[0 0 0],'LineWidth',1.3)

    nexttile
    xlabel('% of cardiac cycle')
    ylabel('E LV')
    hold on
    for i = 1:length(Trange)
        plot(100.*simtimeOrig./T,allHeartbeats.y(i,:,ind.Emax_LV),'color',colors(i,:))
    end
    plot(100*(solTrue.t)./T,solTrue.y(:,ind.Emax_LV),'k-','LineWidth',1.3)
    plot(100*sol.t./sol.t(end),sol.y(:,ind.Emax_LV),'--','Color',[0 0 0],'LineWidth',1.3)

    nexttile
    xlabel('% of cardiac cycle')
    ylabel('E LA')
    hold on
    for i = 1:length(Trange)
        plot(100.*simtimeOrig./T,allHeartbeats.y(i,:,ind.Emax_LA),'color',colors(i,:))
    end
    plot(100*(solTrue.t)./T,solTrue.y(:,ind.Emax_LA),'k-','LineWidth',1.3)
    plot(100*sol.t./sol.t(end),sol.y(:,ind.Emax_LA),'--','Color',[0 0 0],'LineWidth',1.3)

    %in time
    nexttile
    xlabel('time (s)')
    ylabel('MV')
    hold on
    for i = 1:length(Trange)
        plot(allHeartbeats.t(i,:),allHeartbeats.x(i,:,ind.MV),'color',colors(i,:))
    end
    plot(solTrue.t,solTrue.x(:,ind.MV),'k-','LineWidth',1.3)
    plot(sol.t,sol.x(:,ind.MV),'--','Color',[0 0 0],'LineWidth',1.3)

    nexttile
    xlabel('time (s)')
    ylabel('E LV')
    hold on
    for i = 1:length(Trange)
        plot(allHeartbeats.t(i,:),allHeartbeats.y(i,:,ind.Emax_LV),'color',colors(i,:))
    end
    plot(solTrue.t,solTrue.y(:,ind.Emax_LV),'k-','LineWidth',1.3)
    plot(sol.t,sol.y(:,ind.Emax_LV),'--','Color',[0 0 0],'LineWidth',1.3)

    nexttile
    xlabel('time (s)')
    ylabel('E LA')
    hold on
    for i = 1:length(Trange)
        plot(allHeartbeats.t(i,:),allHeartbeats.y(i,:,ind.Emax_LA),'color',colors(i,:))
    end
    plot(solTrue.t,solTrue.y(:,ind.Emax_LA),'k-','LineWidth',1.3)
    plot(sol.t,sol.y(:,ind.Emax_LA),'--','Color',[0 0 0],'LineWidth',1.3)
end



end
function [cost] = costFunction_corr(theta,constants,allparams,simulationOptions,numHeartBeats,data,ind,simulationFunction,dispErrors,doPlot)
% Calculates the cost for the Avatar model. 

%convert the parameters back after optimization
allparams(ind.estParams) = 10.^theta;
theta = allparams;

T = constants(ind.T);

%% SIMULATE THE MODEL FOR THE GIVEN PARAM
simulationOptions.x0 = data.IC;
datatime = data.MV.time';
simtime = sort([datatime,0:0.001:T]);
simtime = unique(simtime);
datatimeinds = find(ismember(simtime,datatime));
try
    simLast = simulate_avatarHEALTH_short(theta,constants,simulationOptions,numHeartBeats,ind,simtime,simulationFunction);
catch myError
    disp("!!Error when simulating!!")
    disp(myError.message)
    cost = 1e99; %Return an extremly large cost if the simulation fails and there is an error.
    return %This command causes the cost function to stop here, and it does not read the subsequent lines
end 


nansum = sum(isnan(simLast.y(:)));
if (simLast.status<0) || (nansum > 0) %if simulation is NaN (ie the simulation failed)
    cost = 1e10;
    doPlot = 0;
    if dispErrors
        disp(['Simulation failed. Nan sum ' num2str(nansum) ', Status ' num2str(simLast.status)])
    end
else
    %% calculate cost for each of the three blood flow measurements
    indsyst = 1:ind.tdiast; %systole
    inddiast = ind.tdiast:length(datatime); %diastole
    
    % which simulations and parameters to be compared with data
    datanames = {'MV','AV','AA','PV'};
    simnames = {'mvCorrVol','avCorrVol','aaCorrVol','pvCorrVol'};%datanames;
    dataparams = {'Emax_LV','Caa','Ctot','Rtot','ELCo'};
    
    costs = zeros(1,length(datanames)+2+length(dataparams));
    residuals = cell(1,length(datanames));
    compSims = cell(1,length(datanames));
    compTime = cell(1,length(datanames));
    compDatas = cell(1,length(datanames));
    compDataUnc = cell(1,length(datanames));
    oscillationPunishment = zeros(1,length(datanames));
    highresSims = cell(1,length(datanames));
    highresTime = cell(1,length(datanames));
    pkLocs = cell(1,length(datanames));

    for i = 1:length(datanames)
        if strcmp('MV',datanames{i})
            dataind = inddiast;
        elseif strcmp('PV',datanames{i}) 
            dataind = 1:length(datatime);
        else
            dataind = indsyst;
            Npeaks = 1;
        end
        simind = datatimeinds(dataind);
        compSims{i} = simLast.y(simind,ind.(simnames{i}))';
        compTime{i} = simLast.t(simind);
        indt = find(ismember(round(simLast.t,2),round(compTime{i},2)));
        highresSims{i} = simLast.y(indt(1):indt(end),ind.(simnames{i}))';
        highresTime{i} = simLast.t(indt(1):indt(end));
        
        compDatas{i} = data.(datanames{i}).mean(dataind);
        compDataUnc{i} =  data.(datanames{i}).sem(dataind);
        residuals{i} = (compSims{i} - compDatas{i}').^2;
        costs(i) = sum(residuals{i}./compDataUnc{i}.^2');
        
        if ~strcmp('MV',datanames{i}) &&  ~strcmp('PV',datanames{i}) %if AV or AC
            [pks,pkLocs{i}] = findpeaks(highresSims{i},'NPeaks',Npeaks+1);
            if length(pks) > Npeaks %if more than expected num of peaks: add cost punishment
                oscillationPunishment(i) = 100;
            end
        end
        
    end
    
    %% Brachial/Aortic pressure
    maxAorticPressure = max(simLast.y(:,ind.brachialPressure));
    minAorticPressure = min(simLast.y(:,ind.brachialPressure));
    
    % diffparam determining difference between brachial and aortic systolic
    % pressure (data from anglo-cardiff trial II
    % http://dx.doi.org/10.1161/HYPERTENSIONAHA.107.105445) 8.6 mean value
    lb = 0.1;ub =18.9;
    diffparam = max(min(data.SBP.mean-maxAorticPressure,ub),lb);%parameter determining the difference between brachial and aortic SBP
    
    estimatedSBP = maxAorticPressure + diffparam;
    estimatedDBP = minAorticPressure;

    costSBP = ((data.SBP.mean - estimatedSBP)^2)/ (data.SBP.sem^2);
    costDBP = ((data.DBP.mean - estimatedDBP)^2)/ (data.DBP.sem^2);
    
    costs(length(datanames)+1) = costSBP;
    costs(length(datanames)+2) = costDBP;

    %% Parameters
    paramcosts = zeros(1,length(data.parameters));
    for p = 1:length(dataparams)
        pval = theta(ind.(dataparams{p}));
        dval = data.parameters.(dataparams{p}).mean;
        dvalsem = data.parameters.(dataparams{p}).sem;
        paramcosts(p) = ((pval - dval)^2) / (dvalsem^2);
    end
    costs(length(datanames)+3:end) = paramcosts;
    
    %% punishment if Emax < Emin (non physiological)
    emaxEminPunishment = 0;
    if theta(ind.Emax_LA) < theta(ind.Emin_LA)
        emaxEminPunishment = 200 + (theta(ind.Emin_LA)-theta(ind.Emax_LA))*10000;
    end
    if theta(ind.Emax_LV) < theta(ind.Emin_LV)
        emaxEminPunishment = emaxEminPunishment + 200 + (theta(ind.Emin_LV)-theta(ind.Emax_LV))*10000;
    end
    
    %% Punishment for offset correction parameters 
    %-40 since the bounds are -40 to 40 but opt 0 to 80.
    offsetPunishment = sum(abs(theta(ind.mvCorr:ind.pvCorr)-40)*20);%increase??
    
    %% Final cost
    cost = sum(costs)+sum(oscillationPunishment) + emaxEminPunishment + offsetPunishment;
end

%% Plotting
if doPlot
    f = figure('Name','Cost plot');
    f.Units = 'Normalized';
    f.OuterPosition = [0 0 1 1 ];
    sgtitle(sprintf('Cost: %0.2f (pure: %0.2f)', cost,sum(costs)))
    hold on
    for i = 1:length(compDatas)
        subplot(2,3,i)
        title(sprintf('%s: %0.2f ',datanames{i},costs(i)))
        hold on
        plot(compTime{i},compDatas{i},'ko-','LineWidth',2)
        errorbar(compTime{i},compDatas{i},compDataUnc{i},'k','LineWidth',1.5)
        plot(highresTime{i},highresSims{i},'r--','LineWidth',2)
        plot(compTime{i},compSims{i},'r*','LineWidth',2)
        for p = 1:length(pkLocs{i})
            plot(highresTime{i}(pkLocs{i}(p)),highresSims{i}(pkLocs{i}(p)),'bo','LineWidth',2)
        end
        legend({'Data','Simulation'})
    end
    
    subplot(2,3,i+1)
    title(sprintf('Params: %0.2f',sum(paramcosts)))
    hold on
    for p = 1:length(dataparams)
        pval = theta(ind.(dataparams{p}));
        dval = data.parameters.(dataparams{p}).mean;
        dsem = data.parameters.(dataparams{p}).sem;
        plot(p,log10(pval),'r*','LineWidth',2)
        errorbar(p,log10(dval),log10(dsem),'ko','LineWidth',2)
    end
    xticks(1:length(dataparams))
   xticklabels(dataparams);
   ylabel('Param value (log10)')
   
   subplot(2,3,i+2)
   title(sprintf('Cost SBP: %0.2f, Cost DBP: %0.2f',costSBP,costSBP))
   hold on
   plot(simLast.t,simLast.y(:,ind.brachialPressure),'r-','LineWidth',2)
   plot(simLast.t(simLast.y(:,ind.brachialPressure) == maxAorticPressure),estimatedSBP,'r*')
   plot(simLast.t(simLast.y(:,ind.brachialPressure) == minAorticPressure),estimatedDBP,'r*')
   yline(data.SBP.mean,'k-','LineWidth',2);
   yline(data.DBP.mean,'k-','LineWidth',2);
   errorbar(simLast.t(simLast.y(:,ind.brachialPressure) == maxAorticPressure),data.SBP.mean,data.SBP.sem,'ko')
   errorbar(simLast.t(simLast.y(:,ind.brachialPressure) == minAorticPressure),data.DBP.mean,data.DBP.sem,'ko')
   xlim([0 simLast.t(end)])
end

end
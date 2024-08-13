function [cost,simLast] = costFunction_chi2(theta,constants,allparams,simulationOptions,data,ind,dispErrors,doPlot)
% Calculates the cost for the Avatar model. 

%convert the parameters back after optimization
allparams(data.meta.paramsToOptimize) = 10.^theta;
theta = allparams;

T = constants(ind.T);


%% SIMULATE THE MODEL FOR THE GIVEN PARAM
simulationOptions.x0 = data.IC;
datatime = data.MV.time';
simtime = sort([datatime,0:0.001:T]);
simtime = unique(simtime);
datatimeinds = find(ismember(simtime,datatime));
doSimulation = str2func(data.meta.simulationFunction);
try
    simLast = doSimulation(theta,constants,simulationOptions,ind,simtime,data.meta.modelName);
catch myError
    disp("!!Error when simulating!!")
    disp(myError.message)
    cost = 1e20; %Return an extremly large cost if the simulation fails and there is an error.
    return %This command causes the cost function to stop here, and it does not read the subsequent lines
end 


nansum = sum(isnan(simLast.x(:)));
if (simLast.status<0) || (nansum > 0) %if simulation is NaN (ie the simulation failed)
    cost = 1e10;
    doPlot = 0;
    if dispErrors
        disp(['Simulation failed. Nan sum ' num2str(nansum) ', Status ' num2str(simLast.status)])
    end
else
    %% calculate cost for each of the blood flow measurements
    indsyst = 1:ind.tdiast; %systole
    inddiast = ind.tdiast:length(datatime); %diastole
    
    % which simulations and parameters to be compared with data
    datanames = {'MV','AV','AA','PV'};
    simnames = {'mvCorrVol','avCorrVol','aaCorrVol','pvCorrVol'};
    dataparams = {'Emax_LV','Caa','Ctot','Rtot','ELCo'};
    
    costs = zeros(1,length(datanames)+2+length(dataparams));
    residuals = cell(1,length(datanames));
    compSims = cell(1,length(datanames));
    compTime = cell(1,length(datanames));
    compDatas = cell(1,length(datanames));
    compDataUnc = cell(1,length(datanames));

    for i = 1:length(datanames)
        if strcmp('MV',datanames{i})
            dataind = inddiast;
        elseif strcmp('PV',datanames{i}) 
            dataind = 1:length(datatime);
        else
            dataind = indsyst;
        end
        simind = datatimeinds(dataind);
        compSims{i} = simLast.y(simind,ind.(simnames{i}))';
        compTime{i} = simLast.t(simind);        
        compDatas{i} = data.(datanames{i}).mean(dataind);
        compDataUnc{i} =  data.(datanames{i}).sem(dataind);
        residuals{i} = compDatas{i}' - compSims{i};
        costs(i) = sum((residuals{i}.^2)./(compDataUnc{i}.^2'));
    end
    
    %% Brachial/Aortic pressure
    maxAorticPressure = max(simLast.y(:,ind.aorticPressure));
    minAorticPressure = min(simLast.y(:,ind.aorticPressure));
    
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
    
    
    %% Final cost
    cost = sum(costs);
end

%% Plotting
if doPlot
    f = figure('Name','Cost plot');
    f.Units = 'Normalized';
    f.OuterPosition = [0 0 1 1 ];
    sgtitle(sprintf('Cost: %0.2f', cost))
    hold on
    for i = 1:length(compDatas)
        subplot(2,3,i)
        title(sprintf('%s: %0.2f',datanames{i},costs(i)))
        hold on
        plot(compTime{i},compDatas{i},'ko-','LineWidth',2)
        errorbar(compTime{i},compDatas{i},compDataUnc{i},'k','LineWidth',1.5)
        plot(compTime{i},compSims{i},'r*','LineWidth',2)
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
   plot(simLast.t,simLast.y(:,ind.aorticPressure),'r-','LineWidth',2)
   plot(simLast.t(simLast.y(:,ind.aorticPressure) == maxAorticPressure),estimatedSBP,'r*')
   plot(simLast.t(simLast.y(:,ind.aorticPressure) == minAorticPressure),estimatedDBP,'r*')
   yline(data.SBP.mean,'k-','LineWidth',2);
   yline(data.DBP.mean,'k-','LineWidth',2);
   errorbar(simLast.t(simLast.y(:,ind.aorticPressure) == maxAorticPressure),data.SBP.mean,data.SBP.sem,'ko')
   errorbar(simLast.t(simLast.y(:,ind.aorticPressure) == minAorticPressure),data.DBP.mean,data.DBP.sem,'ko')
   xlim([0 simLast.t(end)])
end

end
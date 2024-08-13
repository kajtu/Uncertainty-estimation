clear
close all

% Load and format data for 80 subjects, and add the estimated combined errors.
for p = 1:80
    % Load the clinical data
    load(['inputdata_avatarmodel_P' num2str(p) '.mat'])

    % load true data, which contains the combined errors
    load('dataSimulated.mat','dataSimulated')
    trueData = dataSimulated.simRandomSystematic;

    try
        %add mean values & time to the struct containing the final
        %estimation data
        estimationData = struct;
        estimationData.MV.mean = data.MV(:,2);
        estimationData.MV.time = data.MV(:,1);
        estimationData.AV.mean = data.AV(:,2);
        estimationData.AV.time = data.AV(:,1);
        estimationData.AA.mean = data.AC(:,2);
        estimationData.AA.time = data.AC(:,1);
        estimationData.PV.mean = data.PV(:,2);
        estimationData.PV.time = data.PV(:,1);

        estimationData.indtdiast = data.indtdiast;
        estimationData.extra.tdiast = paramdata.tdiast;
        estimationData.extra.indtdiast = paramdata.indtdiast;


        %add estimated combined data uncertainty
        [e.estimatedRRerror] = calcRRsmoothingerror(estimationData); %use empirical function to calculate the RR error
        estimationData.MV.sem = trueData.MV.eRand + trueData.MV.eSyst + abs(e.estimatedRRerror.MV.mean);
        estimationData.AV.sem = trueData.AV.eRand + trueData.AV.eSyst + abs(e.estimatedRRerror.AV.mean);
        estimationData.AA.sem = trueData.AA.eRand + trueData.AA.eSyst + abs(e.estimatedRRerror.AA.mean);
        estimationData.PV.sem = trueData.PV.eRand + trueData.PV.eSyst + abs(e.estimatedRRerror.PV.mean);


        estimationData.SBP.sem = trueData.SBP.sem;%inherit the sem
        estimationData.DBP.sem = trueData.DBP.sem;%inherit the sem
        estimationData.extra.ESV.sem = trueData.extra.ESV.eRand;

        estimationData.SBP.mean = paramdata.SBP;
        estimationData.DBP.mean = paramdata.DBP;
        estimationData.extra.ESV.mean = paramdata.ESV_seg;

        estimationData.extra.EDV = NaN;
        estimationData.extra.LaESV = NaN;

        % information for simulating the data
        estimationData.meta.simulationFunction = trueData.meta.simulationFunction;
        estimationData.meta.modelName = trueData.meta.modelName;
        estimationData.meta.paramsToOptimize = trueData.meta.paramsToOptimize;

        % Add sem for all data-based parameters
        estimationData.extra.T.mean = paramdata.T;
        estimationData.parameters.Caa.mean = paramdata.Caa;
        estimationData.parameters.ELCo.mean = paramdata.ElCo;

        [estimationData.parameters.Ctot.mean,estimationData.parameters.Rtot.mean,estimationData.parameters.Emax_LV.mean] = calculateDatabasedParameters(estimationData);

        pnames = fieldnames(trueData.parameters);
        for i = 1:length(pnames)
            estimationData.parameters.(pnames{i}).sem = trueData.parameters.(pnames{i}).sem;%inherit the sem
        end

        %% Optional: create plot of the resulting estimation data to make sure it is reasonable
        figure()
        nexttile
        title('mv')
        errorbar(estimationData.MV.time,estimationData.MV.mean,estimationData.MV.sem)
        nexttile
        title('av')
        errorbar(estimationData.AV.time,estimationData.AV.mean,estimationData.AV.sem)
        nexttile
        title('aa')
        errorbar(estimationData.AA.time,estimationData.AA.mean,estimationData.AA.sem)
        nexttile
        title('pv')
        errorbar(estimationData.PV.time,estimationData.PV.mean,estimationData.PV.sem)
        nexttile
        title('parameters')
        hold on
        for i = 1:length(pnames)
            errorbar(i,estimationData.parameters.(pnames{i}).mean,estimationData.parameters.(pnames{i}).sem,'*');
        end

        nexttile
        hold on
        errorbar(1,estimationData.SBP.mean,estimationData.SBP.sem)
        errorbar(2,estimationData.DBP.mean,estimationData.DBP.sem)
        title('bp')

        %% Save the created estimation data for this subject
        save(['dataP' num2str(p) '.mat'],'estimationData')

    catch
        disp(['couldnt do p ' num2str(p)])
    end

end

%% Plot selected data
load('dataP78.mat','estimationData')
cdata = estimationData;
load('dataP33.mat','estimationData')
pdata = estimationData;

figure()
nexttile
title('mv')
hold on
errorbar(cdata.MV.time./cdata.MV.time(end),cdata.MV.mean,cdata.MV.sem)
errorbar(pdata.MV.time./pdata.MV.time(end),pdata.MV.mean,pdata.MV.sem)

nexttile
title('av')
hold on
errorbar(cdata.AV.time./cdata.MV.time(end),cdata.AV.mean,cdata.AV.sem)
errorbar(pdata.AV.time./pdata.MV.time(end),pdata.AV.mean,pdata.AV.sem)

nexttile
hold on
title('aa')
errorbar(cdata.AA.time./cdata.MV.time(end),cdata.AA.mean,cdata.AA.sem)
errorbar(pdata.AA.time./pdata.MV.time(end),pdata.AA.mean,pdata.AA.sem)

nexttile
hold on
title('pv')
errorbar(cdata.PV.time./cdata.MV.time(end),cdata.PV.mean,cdata.PV.sem)
errorbar(pdata.PV.time./pdata.MV.time(end),pdata.PV.mean,pdata.PV.sem)

nexttile
title('parameters')
hold on
for i = 1:length(pnames)
    errorbar(i,cdata.parameters.(pnames{i}).mean,cdata.parameters.(pnames{i}).sem,'*');
    errorbar(i+0.2,pdata.parameters.(pnames{i}).mean,pdata.parameters.(pnames{i}).sem,'*');
end
xticks(1:length(pnames))
xticklabels(pnames)

nexttile
hold on
errorbar(1,cdata.SBP.mean,cdata.SBP.sem)
errorbar(2,cdata.DBP.mean,cdata.DBP.sem)
errorbar(1.2,pdata.SBP.mean,pdata.SBP.sem)
errorbar(2.2,pdata.DBP.mean,pdata.DBP.sem)
title('bp')

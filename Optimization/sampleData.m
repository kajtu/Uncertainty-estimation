function [estimationData,e] = sampleData(trueData,e,n)
% Function to take the trueData and adding
% errors defined in the struct e. If e is not provided,
% e is defined and samples from the error distribution defined by the
% trueData struct.

% trueData is a struct containing the "true" data including time,
% meanvalues, and the distribution of the measurement uncertainty (sem)
% estimationData inherits the sigma values from trueData

estimationData = struct;

%create the random errors, otherwise use the one provided
if nargin<2 || isempty(e)
    % RR variation: average of 400 heartbeats
    rrNumber = trueData.meta.RRnum;
    if rrNumber > 1
        if nargin>2 && ~isempty(n)
            [~,eloaded] = loadSampledMeasurementError(n,'simRandomSystematic');
            disp('loaded RR smoothing')
            e.RR.MV = eloaded.RR.MV;
            e.RR.AV = eloaded.RR.AV;
            e.RR.AA = eloaded.RR.AA;
            e.RR.PV = eloaded.RR.PV;

            load('dataSimulated.mat','simulatedDataTrue');
            sind = 1:length(e.RR.MV);
            sol.t =  trueData.MV.time;
            sol.x(:,simulatedDataTrue.ind.MV) = e.RR.MV + trueData.MV.mean;
            sol.x(:,simulatedDataTrue.ind.AV) = e.RR.AV + trueData.AV.mean;
            sol.x(:,simulatedDataTrue.ind.AA) = e.RR.AA + trueData.AA.mean;
            sol.x(:,simulatedDataTrue.ind.PV) = e.RR.PV + trueData.PV.mean;
        else
            doPlot = 0;
            randomDistribution = 1;
            load('dataSimulated.mat','simulatedDataTrue');
            fprintf('Simulating %d heartbeats...\n',rrNumber)
            [sol,solTrue,~] = simulate_avatar_RR(simulatedDataTrue.allParameters,simulatedDataTrue.constants,simulatedDataTrue.options,simulatedDataTrue.ind,simulatedDataTrue.simtime,simulatedDataTrue.modelName,rrNumber,randomDistribution,doPlot);
            sind = find(ismember(round(solTrue.t,4),round(trueData.MV.time,4)));
            [~,sindunique] = unique(solTrue.t(sind));
            sind = sind(sindunique);

            e.RR.MV = sol.x(sind,simulatedDataTrue.ind.MV) - trueData.MV.mean;
            e.RR.AV = sol.x(sind,simulatedDataTrue.ind.AV) - trueData.AV.mean;
            e.RR.AA = sol.x(sind,simulatedDataTrue.ind.AA) - trueData.AA.mean;
            e.RR.PV = sol.x(sind,simulatedDataTrue.ind.PV) - trueData.PV.mean;
        end
    else
        e.RR.MV = 0;
        e.RR.AV = 0;
        e.RR.AA = 0;
        e.RR.PV = 0;
    end

    % RR variation: add the estimated size of the RR variation error to the
    % sigma/sem based on the simulated mean curve
    if trueData.meta.addRRtoSem
        rrdata.MV.mean = sol.x(sind,simulatedDataTrue.ind.MV);
        rrdata.AV.mean = sol.x(sind,simulatedDataTrue.ind.AV);
        rrdata.AA.mean = sol.x(sind,simulatedDataTrue.ind.AA);
        rrdata.PV.mean = sol.x(sind,simulatedDataTrue.ind.PV);
        rrdata.MV.time = sol.t(sind);
        rrdata.AV.time = sol.t(sind);
        rrdata.AA.time = sol.t(sind);
        rrdata.PV.time = sol.t(sind);
        rrdata.indtdiast = trueData.extra.indtdiast;
        warning('off')
        [e.estimatedRRerror] = calcRRsmoothingerror(rrdata); 
        warning('on')
    end

    % sample timepoint-specific random error for the blood flow curves
    e.mv.random = normrnd(0,trueData.MV.eRand);
    e.av.random = normrnd(0,trueData.AV.eRand);
    e.aa.random = normrnd(0,trueData.AA.eRand);
    e.pv.random = normrnd(0,trueData.PV.eRand);

    % sample systematic error for all timepoints for the blood flow curves
    e.mv.systematic = normrnd(0,trueData.MV.eSyst);
    e.av.systematic = normrnd(0,trueData.AV.eSyst);
    e.aa.systematic = normrnd(0,trueData.AA.eSyst);
    e.pv.systematic = normrnd(0,trueData.PV.eSyst);

    % sample random and systematic error for the parameters
    e.SBP.random  = normrnd(0,trueData.SBP.eRand);%mu=0, sigma=.eRand
    e.DBP.random  = normrnd(0,trueData.DBP.eRand);
    e.SBP.systematic = normrnd(0,trueData.SBP.eSyst);
    e.DBP.systematic = normrnd(0,trueData.DBP.eSyst);

    pnames = fieldnames(trueData.parameters);
    for i = 1:length(pnames)
        e.(pnames{i}).random = normrnd(0,trueData.parameters.(pnames{i}).eRand);
        e.(pnames{i}).systematic = normrnd(0,trueData.parameters.(pnames{i}).eSyst);
    end

    % sample error for ventricular volume (needed for Emax_LV)
    e.ESV.random = normrnd(0,trueData.extra.ESV.eRand);
    e.ESV.systematic = normrnd(0,trueData.extra.ESV.eSyst);
end

% add both systematic and random errors to the blood flow curves
estimationData.MV.mean = trueData.MV.mean + e.mv.random + e.mv.systematic + e.RR.MV;
estimationData.AV.mean = trueData.AV.mean + e.av.random + e.av.systematic + e.RR.AV;
estimationData.AA.mean = trueData.AA.mean + e.aa.random + e.aa.systematic + e.RR.AA;
estimationData.PV.mean = trueData.PV.mean + e.pv.random + e.pv.systematic + e.RR.PV;

%inherit the sem
if trueData.meta.addRRtoSem %if RR variation should be added to the sem, add it (different for each sample)
    estimationData.MV.sem = trueData.MV.eRand + trueData.MV.eSyst + abs(e.estimatedRRerror.MV.mean);
    estimationData.AV.sem = trueData.AV.eRand + trueData.AV.eSyst + abs(e.estimatedRRerror.AV.mean);
    estimationData.AA.sem = trueData.AA.eRand + trueData.AA.eSyst + abs(e.estimatedRRerror.AA.mean);
    estimationData.PV.sem = trueData.PV.eRand + trueData.PV.eSyst + abs(e.estimatedRRerror.PV.mean);
else %otherwise, inherit it
    estimationData.MV.sem = trueData.MV.sem;
    estimationData.AV.sem = trueData.AV.sem;
    estimationData.AA.sem = trueData.AA.sem;
    estimationData.PV.sem = trueData.PV.sem;
end

%inherit the time
estimationData.MV.time = trueData.MV.time;
estimationData.AV.time = trueData.AV.time;
estimationData.AA.time = trueData.AA.time;
estimationData.PV.time = trueData.PV.time;

%% add random and systematic error for the brachial pressure
estimationData.SBP.mean = trueData.SBP.mean + e.SBP.random + e.SBP.systematic;
estimationData.DBP.mean = trueData.DBP.mean + e.DBP.random + e.DBP.systematic;
estimationData.SBP.sem = trueData.SBP.sem;%inherit the sem
estimationData.DBP.sem = trueData.DBP.sem;%inherit the sem

%% add error for the ventricular segmentation volume
estimationData.extra.ESV.mean  = trueData.extra.ESV.mean + e.ESV.random + e.ESV.systematic; % needed for EmaxLV
estimationData.extra.ESV.sem = e.ESV.random;

%% add extra information needed for simulation and comparison with the simulation
estimationData.extra.indtdiast = trueData.extra.indtdiast; %used to only fit to data during systole/diastole (TODO: calculate for each sampled data?)
estimationData.extra.tdiast = trueData.extra.tdiast;
estimationData.extra.T.mean = trueData.extra.T.mean;
estimationData.extra.LaESV = NaN; %needed for initial conditions, but set to litterature value if NaN
estimationData.extra.EDV = NaN; %needed for initial conditions, but set to litterature value if NaN

% information for simulating the data
estimationData.meta.simulationFunction = trueData.meta.simulationFunction;
estimationData.meta.modelName = trueData.meta.modelName;
estimationData.meta.paramsToOptimize = trueData.meta.paramsToOptimize;

%% add random and systematic error for the parameters
% ELCo and Caa mean values are sampled
pnames = {'ELCo','Caa'};
for i = 1:length(pnames)
    estimationData.parameters.(pnames{i}).mean = trueData.parameters.(pnames{i}).mean + e.(pnames{i}).random + e.(pnames{i}).systematic;%add the error
end

% Ctot, Rtot, and Emax_LV mean values are calculated from all other sampled data
[estimationData.parameters.Ctot.mean,estimationData.parameters.Rtot.mean,estimationData.parameters.Emax_LV.mean] = calculateDatabasedParameters(estimationData);

% Add sem for all data-based parameters
pnames = fieldnames(trueData.parameters);
for i = 1:length(pnames)
    estimationData.parameters.(pnames{i}).sem = trueData.parameters.(pnames{i}).sem;%inherit the sem
end

end
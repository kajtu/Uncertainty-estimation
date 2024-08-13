%% createTrueData_simulated
clear
% Script to create estimation data based on a model simulation as "true" data
% also add the distribution of errors, that can be sampled using
% createSampledMeasurementError (in sampleData)

% NOTE: data names cannot contain underscore!

% add paths
basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');

%add all project folders to the matlab path
addpath(genpath(basefolder))

% setup AMICI toolbox
run(fullfile(basefolder, 'Requirements', 'AMICI-0.10.11_SS_eventFix', 'matlab', 'installAMICI.m'))


%% Create the "true" simulation data

%setup
% Simulation settings
% sensi: sensitivity order - 0 = no calculation of sensitivities (default)
% maxsteps: maximum number of integration steps
simulationoptions = amioption('sensi',0,'maxsteps',1e4);
% set simulation tolerances
simulationoptions.atol = 1e-16;
simulationoptions.rtol = 1e-8;

load('trueParameters.mat','trueParams','trueConstants','truedata','ind','paramdata','extradata')
paramNamesOld = {'Cpvc' 'Rpu' 'Rpv' 'Lpv' 'Rtot' 'Ctot' 'ElCo' 'Caa' 'Emax_LA' 'Emax_LV' 'Emin_LA' 'Emin_LV' 'Lao' 'Lav'...
    'Lmv' 'Ppu' 'Rao' 'Rmv' 'k_diast_LA' 'k_diast_LV' 'k_syst_LA'...
    'k_syst_LV' 'm1_LA' 'm1_LV' 'm2_LA' 'm2_LV' 'onset_LA' 'onset_LV'};
constantsNamesOld = {'aaCorr' 'avCorr' 'mvCorr' 'tdiast' 'Ks_LA' 'Ks_LV'...
    'V0_LA' 'V0_LV' 'RLAvisc' 'RLVvisc' 'Raa' 'Rpc' 'Rpvc' 'T' 'rho_blood' 'norm_factor_LA' 'norm_factor_LV'};

%load parameter names and indexes
T = trueConstants(ind.T);
simulatedDataTrue.extra.T.mean = T;
simulatedDataTrue.extra.indtdiast = paramdata.indtdiast;
simulatedDataTrue.extra.tdiast = paramdata.tdiast;
simulatedDataTrue.parameters.ELCo.mean = trueParams(ind.ElCo);
simulatedDataTrue.parameters.Caa.mean = trueParams(ind.Caa);
simulatedDataTrue.parameters.Ctot.mean = trueParams(ind.Ctot);
simulatedDataTrue.parameters.Rtot.mean = trueParams(ind.Rtot);
simulatedDataTrue.parameters.Emax_LV.mean = trueParams(ind.Emax_LV);
oldind = ind;
[~,~,paramNames,constantsNames,units,ind] = loadParameters(simulatedDataTrue);

%convert the true params to the new order or params and constants
trueConstants = trueConstants(4:end); %remove corr constants
trueParams = [trueParams,40,40,40,40]; % add correction parameters

%simulate
simulationoptions.x0 = truedata.IC;
datatime = truedata.time;
simtime = sort([datatime,0:0.001:T]);
simtime = unique(simtime);
datatimeinds = find(ismember(simtime,datatime));
cd ../Simulation
warning('off')
[sol] = simulate_avatar_ss(trueParams,trueConstants,simulationoptions,ind,simtime,'avatar_corr');
warning('on')
cd ../Data
%set the simulated data based on the simulation
simulatedDataTrue.MV.time = sol.t;
simulatedDataTrue.AV.time = sol.t;
simulatedDataTrue.AA.time = sol.t;
simulatedDataTrue.PV.time = sol.t;
simulatedDataTrue.MV.mean = sol.x(:,ind.MV);
simulatedDataTrue.AV.mean = sol.x(:,ind.AV);
simulatedDataTrue.AA.mean = sol.x(:,ind.AA);
simulatedDataTrue.PV.mean = sol.x(:,ind.PV);
simulatedDataTrue.extra.SV_MV.mean = trapz(sol.t,sol.x(:,ind.MV));
simulatedDataTrue.extra.SV_AV.mean = trapz(sol.t,sol.x(:,ind.AV));
simulatedDataTrue.extra.SV_AA.mean = trapz(sol.t,sol.x(:,ind.AA));
simulatedDataTrue.extra.EDV.mean = max(sol.x(:,ind.LV));
simulatedDataTrue.extra.ESV.mean = min(sol.x(:,ind.LV));
simulatedDataTrue.extra.LaESV = max(sol.y(:,ind.Vla));
simulatedDataTrue.Vlv = sol.x(:,ind.LV);
d.SBP.mean = truedata.SBP;
[simulatedDataTrue.SBP.mean,simulatedDataTrue.DBP.mean] = brachialpressure(sol,ind,d);
simulatedDataTrue.SBP.time = NaN;simulatedDataTrue.DBP.time = NaN;

%% Calculate RR smoothing
rrNumber = 400;
randomDistribution = 0;
fprintf('Simulating %d heartbeats...\n',rrNumber)
[sol,solTrue,~] = simulate_avatar_RR(trueParams,trueConstants,simulationoptions,ind,simtime,'avatar_corr',rrNumber,randomDistribution,0); 
rrdata.MV.mean = sol.x(:,ind.MV);
rrdata.AV.mean = sol.x(:,ind.AV);
rrdata.AA.mean = sol.x(:,ind.AA);
rrdata.PV.mean = sol.x(:,ind.PV);
rrdata.MV.time = sol.t;
rrdata.AV.time = sol.t;
rrdata.AA.time = sol.t;
rrdata.PV.time = sol.t;
[~,rrdata.indtdiast] = min(abs(sol.t-paramdata.tdiast));
[estimatedRRerror] = calcRRsmoothingerror(rrdata); %re-estimated for each sampled data in SampleData (createSampledMeasurementError)

%% Add uncertainty: only random error included in sigma (sem)
exp = simulatedDataTrue;
exp.meta.timeResolution = 40;
exp = setTemporalResolution(exp,exp.meta.timeResolution);

% Blood flow curves: standard deviation of the random error
exp.MV.eRand = abs(exp.MV.mean).*0.0539; % 5.4% of flow  - different in each timepoint
exp.MV.eSyst= 0.9*7; % 0.9 cm/s * vessel area  (doi:10.1016/j.jsha.2010.07.007, https://doi.org/10.1152/ajpheart.00269.2004 ) - same for all timepoints
exp.AV.eRand = abs(exp.AV.mean).*0.0539; % 5.4% of flow
exp.AV.eSyst= 0.9*5; % 0.9 cm/s * vessel area (doi: 10.4137/CMC.S15716, doi:10.1053/euhj.2001.2782,  pi*(3/2)^2)
exp.AA.eRand = abs(exp.AA.mean).*0.0539; % 5.4% of flow
exp.AA.eSyst= 0.9*7;% 0.9 cm/s * vessel area (doi:10.1053/euhj.2001.2782   pi*(2.5/2)^2)
exp.PV.eRand = abs(exp.PV.mean).*0.0806; % 12% of flow
exp.PV.eSyst = 0.9*5.3;% 0.9 cm/s * vessel area (4*pi*(1.3/2)^2 = 5.3 cm^2  (Wittkampf et al., 2003; Kim et al., 2005))

% where the flow is close to 0, the % error is too small. Correct this:
exp.MV.eRand(exp.MV.eRand < mean(exp.MV.eRand)) = mean(exp.MV.eRand);
exp.AV.eRand(exp.AV.eRand < mean(exp.AV.eRand)) = mean(exp.AV.eRand);
exp.AA.eRand(exp.AA.eRand < mean(exp.AA.eRand)) = mean(exp.AA.eRand);
exp.PV.eRand(exp.PV.eRand < mean(exp.PV.eRand)) = mean(exp.PV.eRand); 

% add the error to sigma (only the random error)
exp.MV.sem = exp.MV.eRand; 
exp.AV.sem = exp.AV.eRand;
exp.AA.sem = exp.AA.eRand;
exp.PV.sem = exp.PV.eRand;

% Data-based parameters
exp.parameters.Caa.eRand = simulatedDataTrue.parameters.Caa.mean.*0.15; % 15% error
exp.parameters.ELCo.eRand = simulatedDataTrue.parameters.ELCo.mean.*0.12; % 12% error
exp.parameters.Ctot.eRand = simulatedDataTrue.parameters.Ctot.mean.*0.20; % 25%, calculated from other variables
exp.parameters.Rtot.eRand = simulatedDataTrue.parameters.Rtot.mean.*0.082; % 12%, calculated from other variables
exp.parameters.Emax_LV.eRand = simulatedDataTrue.parameters.Emax_LV.mean.*0.077; % 10%, calculated from other variables

exp.parameters.Caa.eSyst = 0; % no systematic errors
exp.parameters.ELCo.eSyst = 0; % no systematic errors
exp.parameters.Ctot.eSyst = 0; % no systematic errors added, the parameter is calculated from other measurements
exp.parameters.Rtot.eSyst = 0; % no systematic errors added, the parameter is calculated from other measurements
exp.parameters.Emax_LV.eSyst = 0; % no systematic errors added, the parameter is calculated from other measurements

exp.parameters.Caa.sem = exp.parameters.Caa.eRand;
exp.parameters.ELCo.sem = exp.parameters.ELCo.eRand;
exp.parameters.Ctot.sem = exp.parameters.Ctot.eRand;
exp.parameters.Rtot.sem = exp.parameters.Rtot.eRand;
exp.parameters.Emax_LV.sem = exp.parameters.Emax_LV.eRand;

% Brachial pressure
exp.SBP.eRand = 8; % 8 mmHg error
exp.SBP.eSyst= 0;
exp.SBP.sem = exp.SBP.eRand;

exp.DBP.eRand = 8; % 8 mmHg error
exp.DBP.eSyst = 0;
exp.DBP.sem = exp.DBP.eRand;

% LV segmentation
exp.extra.ESV.eRand = simulatedDataTrue.extra.ESV.mean * 0.05; % 5% random error
exp.extra.ESV.eSyst= 0;
exp.extra.ESV.sem = exp.extra.ESV.eRand;

% Other information
exp.meta.RRnum = 400; %average of 400 random heartbeats
exp.meta.addRRtoSem = 0; %do not add the systematic RR error to sigme
exp.meta.description = 'only random errors in sigma';
exp.meta.paramsToOptimize = 1:length(paramNames);%include the corr parameters
exp.meta.simulationFunction = 'simulate_avatar_correction';
exp.meta.modelName = 'avatar_corr';
expname = 'simRandomOnly';
dataSimulated.(expname) = exp;

%% Add uncertainty: both random and systematic errors in the sigma/sem
%calculate rr error
sind = find(ismember(solTrue.t,exp.MV.time));
[~,sindunique] = unique(solTrue.t(sind));
sind = sind(sindunique);
datanames = {'MV','AV','AA','PV'};
for n = 1:length(datanames)
     estimatedRRerror.(datanames{n}).mean = estimatedRRerror.(datanames{n}).mean(sind)';
end

%all errors are the same as above, it is just the sigma and the simulation
%function that changes (plus which parameters that are estimated)
exp.MV.sem = exp.MV.eRand + exp.MV.eSyst + abs(estimatedRRerror.MV.mean');
exp.AV.sem = exp.AV.eRand + exp.AV.eSyst + abs(estimatedRRerror.AV.mean');
exp.AA.sem = exp.AA.eRand + exp.AA.eSyst + abs(estimatedRRerror.AA.mean');
exp.PV.sem = exp.PV.eRand + exp.PV.eSyst + abs(estimatedRRerror.PV.mean');

exp.meta.simulationFunction = 'simulate_avatar_ss';
exp.meta.modelName = 'avatar_corr';
exp.meta.paramsToOptimize = 1:find(strcmp(paramNames,'onset_LV'));%do not include the corr parameters
exp.meta.addRRtoSem = 1; %DO add the systematic RR error to sigma
expname = 'simRandomSystematic';
dataSimulated.(expname) = exp;


%% Save the full data structure
% add extra information
simulatedDataTrue.allSimulations = sol;
simulatedDataTrue.allParameters = trueParams;
simulatedDataTrue.constants = trueConstants;
simulatedDataTrue.ind = ind;
simulatedDataTrue.options = simulationoptions;
simulatedDataTrue.simtime = simtime;
simulatedDataTrue.simFunc = 'simulate_avatar_corr';
simulatedDataTrue.modelName = 'avatar_corr';

save('dataSimulated.mat','dataSimulated','simulatedDataTrue')


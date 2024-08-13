%runParamEstimationSeveral_ESS_asNSC
% This scripts runs parameter estimation in several iterations, similar to
% what was done in the manuscript. However, in the manuscript the
% parallelisation was done directly in bash scripts that were run on the
% Swedish National Supercomputing Center (NSC). 
clear
close all
    
rng('shuffle')

%% Set variables
experimentNames = 'simRandomOnly';%simRandomOnly or  simRandomSystematic
experimentRange = '1';
numExperiments = 1;
startN = 0;
numberOfSamples=1;
method = 'ESS';
continueMCMC = 0;
doParamEst=1;
datename = datestr(now,'yymmdd-HHMMSS');
addpath ../Simulation
% Setup(experimentNames,doParamEst,datename,experimentRange,method,continueMCMC,numExperiments,startN)
restoredefaultpath

dataName = 'dataP78';

%% Run optimizations
tic
randomstart = 1; 
disp('***Starting optimization run 1 - random starts***')
for i = 1:(40*numExperiments) %40 per experiment
% parfor i = 1:4%(40*numExperiments) %40 per experiment
    % EstimateParametersESS(i,datename,experimentNames,experimentRange,randomstart,numberOfSamples,startN) %for sampled data
     EstimateParametersESS_realdata(i,datename,dataName,randomstart) %for real data
end
randomstart = 0; 
disp('***Starting optimization run 2 - from previous results***')
for i = 1:(10*numExperiments) %10 per experiment, startguess from previous results
    % EstimateParametersESS(i,datename,experimentNames,experimentRange,randomstart,numberOfSamples,startN) %for sampled data
    EstimateParametersESS_realdata(i,datename,dataName,randomstart) % for real data
end

%% Save information about the time of optimization
timetotalOpt = toc
timetotalOpt_minutes = timetotalOpt/60
timetotalOpt_hours = timetotalOpt_minutes/60

folder= ['..' filesep 'Parameters' filesep 'ESS' filesep];
save([folder 'timeOpt_' method '_' datename '.mat'],'timetotalOpt','timetotalOpt_minutes','timetotalOpt_hours','experimentNames','experimentRange')




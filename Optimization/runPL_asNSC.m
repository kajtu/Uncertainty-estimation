%runPL_asNSC
% This scripts runs profile likelihood in several iterations, similar to
% what was done in the manuscript. However, in the manuscript the
% parallelisation was done directly in a bash scripts that were run on the
% Swedish National Supercomputing Center (NSC). 

clear
close all
    
rng('shuffle')

%% Set variables
experimentNames = 'simRandomOnly';%simRandomOnly or simRandomSystematic
numExperiments = 100;
startN = 0;
% numberOfSamples=numExperiments;
method = 'PL';
paramToEstimate = '';
doParamEst=1;
continueFromLimit=1;
datename = datestr(now,'yymmdd-HHMMSS');
addpath ../Simulation
% Setup(experimentNames,doParamEst,datename,experimentRange,method,continueMCMC,numExperiments,startN)
restoredefaultpath

dataName = 'dataP33';
numberOfSamples = 32;
addRandomness = 0;

%% Run optimizations
disp('***Starting optimization Profile Likelihood ***')
for i = numberOfSamples % choose number of iterations here
% parfor i = 1:4 % use a parfor loop if parallelising the runs
    % For real data, classic method:
    EstimatePL_realdata_classic(i,datename,dataName,continueFromLimit,paramToEstimate,0)

    % Call any PL script here, such as, for sample data:
    % EstimatePL(i,datename,experimentNames,experimentRange,numberOfSamples,startN,paramToEstimate) 
    % EstimatePL_refine(i,datename,experimentNames,experimentRange,numberOfSamples,startN,paramToEstimate)
    % EstimatePL_smooth(i,datename,experimentNames,experimentRange,numberOfSamples,startN,paramToEstimate)

    % For real data and pre-defined parameter steps:
    % EstimatePL_realdata(i,datename,dataName,numberOfSamples,paramToEstimate,1,0,0,0,0,1)
    % EstimatePL_realdata(i,datename,dataName,numberOfSamples,paramToEstimate,0,1,16,16,0,0)
    % Argument names: loopNum,dateStart,dataName,numSamples,paramToEstimate,addRandomness,increasebound,lowersteps,uppersteps,refine,smooth
end


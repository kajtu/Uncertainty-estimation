%% collectPLparams for several subjects and parameters.
% For each subject and biomarker, all PL steps together with the best cost for each step is saved in the file
% allcostsLoaded_xx.mat in the subject folder.

% This script should be run before starting a new round of PL, since the PL
% script loads the allcostsLoaded_xx.mat file and not the individual
% folders.

%How to run at interactive node:
% interactive -n1 -t 00:50:00 or  interactive -n1 -t 00:55:00 --reservation=devel
% cd '/proj/cardiovascular_modeling/users/x_kajtu/Uncertainty-estimation/Optimization'
% module load MATLAB/2023b-bdist
% MATLAB_0='matlab -nodesktop -nodisplay -singleCompThread'
% ${MATLAB_0} -r "collectAllPLparams; exit"

%% Setup paths
basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');
% % add all project folders to the matlab path
addpath(genpath(fullfile(basefolder,'Optimization')))
addpath(genpath(fullfile(basefolder,'Data')))
addpath(genpath(fullfile(basefolder,'Simulation')))
addpath(genpath(fullfile(basefolder,'Tools')))
% setup AMICI and MEIGO toolboxes
run(fullfile(basefolder, 'Requirements', 'AMICI-0.10.11_SS_eventFix', 'matlab', 'installAMICI.m'))

%% Setup variable names
paramNamesPlot = {'Caa', 'Rao','m2_LV','Cpvc','Emax_LA','k_syst_LV'}; % define which biomarkers

experimentNames = {'E_dataP78','E_dataP1','E_dataP3','E_dataP24','E_dataP36','E_dataP33'}; %define which datasets

paramNames = {'Cpvc' 'Rpu' 'Rpv' 'Lpv' 'Rtot' 'Ctot' 'ELCo' 'Caa' 'Emax_LA' 'Emax_LV' 'Emin_LA' 'Emin_LV' 'Lao' 'Lav'...
    'Lmv' 'Ppu' 'Rao' 'Rmv' 'k_diast_LA' 'k_diast_LV' 'k_syst_LA'...
    'k_syst_LV' 'm1_LA' 'm1_LV' 'm2_LA' 'm2_LV' 'onset_LA' 'onset_LV'};
resultsFolder = fullfile(basefolder,'Parameters');

%% Run collectPLparams
disp('Setup done')

for e = 1:length(experimentNames)
    experiment = experimentNames{e}(3:end);
    disp(experiment)
    collectPLparams(fullfile(resultsFolder,'PL',['E_' experiment]),paramNamesPlot,paramNames,fullfile(resultsFolder,'ESS',['E_' experiment '*']));
end

disp('collectAllPLparams done!')
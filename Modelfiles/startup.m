%% Adds paths to toolboxes needed for simulations
basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');


addpath(fullfile(basefolder,'Requirements','PESTO-1.1.0'))
addpath(fullfile(basefolder,'Requirements','AMICI-0.10.11_SS_eventFix','matlab'))
addpath(fullfile(basefolder,'Requirements','MEIGO'))

run(fullfile(basefolder,'Requirements','AMICI-0.10.11_SS_eventFix','matlab','installAMICI.m'))
run(fullfile(basefolder,'Requirements','MEIGO','install_MEIGO.m'))

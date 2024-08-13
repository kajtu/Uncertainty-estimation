function []=Setup(experimentNames,doParamEst,datename,experimentRange,method,continueMCMC,numberOfSamples,startN)

basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');

%add all project folders to the matlab path
addpath(genpath(fullfile(basefolder,'Optimization')))
addpath(genpath(fullfile(basefolder,'Data')))
addpath(genpath(fullfile(basefolder,'Modelfiles')))
addpath(genpath(fullfile(basefolder,'Requirements')))
addpath(genpath(fullfile(basefolder,'Simulation')))

% setup AMICI and MEIGO toolboxes
run(fullfile(basefolder, 'Requirements', 'AMICI-0.10.11_SS_eventFix', 'matlab', 'installAMICI.m'))
run(fullfile(basefolder, 'Requirements', 'MEIGO', 'install_MEIGO.m'))

%METHOD:
% inverse profile likelihood: 'invPL'
% inverse prediction profile likelihood: 'PPLinv'
% profile likelihood: 'PL'
% "normal" parameter optimization: 'opt'


if doParamEst
    %% Create necessary folders
    if ~exist('Log','dir') % You "need" to have a log folder
        mkdir('Log')
    end
    

    if ischar(experimentNames) || isstring(experimentNames)
        experimentNames = {experimentNames};
    end
    if ~strcmp(experimentNames{1}(1:3),'sim')% if not simulated data
        % find the wanted experiments
        load(fullfile(basefolder,'Data','data.mat'),'data')
        [~,experimentNames] = extractData(experimentNames,experimentRange,data);
    end

    %% create folders to save results in
    % Find names of results folder depending on what type of param
    % estimation method
    if nargin > 5
        if strcmp(method,'PPLinv')
            resultsfolder = ['Parameters' filesep 'PPLinv'];
        elseif strcmp(method,'PPL')
            resultsfolder = ['Parameters' filesep 'PPL'];
        elseif strcmp(method,'PL')
            resultsfolder = ['Parameters' filesep 'PL'];
        elseif strcmp(method,'PLinv')
            resultsfolder = ['Parameters' filesep 'PLinv'];
        elseif strcmp(method,'MCMC')
            resultsfolder = ['Parameters' filesep 'MCMC'];
        elseif strcmp(method,'ESS')
            resultsfolder = ['Parameters' filesep 'ESS'];
        else
            resultsfolder = 'Parameters';
        end
    else
        resultsfolder = 'Parameters';
    end
    
    
    disp(['Setup: Number of experiments: ' num2str(length(experimentNames))])
    if continueMCMC
        disp('continuing on old mcmc run')
    else

        if strcmp(experimentNames{1}(1:3),'sim')
            for n = 1:numberOfSamples
                num = n+startN;
                if strcmp(method,'PPLinv') || strcmp(method,'PLinv') || strcmp(method,'PL') || strcmp(method,'PPL')
                    folder = fullfile(basefolder, resultsfolder, ['E_' experimentNames{1} num2str(num)]);
                else
                    folder = fullfile(basefolder, resultsfolder, ['E_' experimentNames{1} num2str(num) '_' datename]);
                end
                if isempty(dir(folder))
                    mkdir(folder)
                    disp(['Setup: created folder ' folder])
                else
                    disp(['Setup: using existing folder: ' folder])
                end
            end
        else
            for e = 1:length(experimentNames)
                if strcmp(method,'PPLinv') || strcmp(method,'PLinv')  || strcmp(method,'PL') || strcmp(method,'PPL')
                    folder = fullfile(basefolder, resultsfolder, ['E_' experimentNames{e}]);
                else
                    folder = fullfile(basefolder, resultsfolder, ['E_' experimentNames{e} '_' datename]);
                end
                if isempty(dir(folder))
                    mkdir(folder)
                    disp(['Setup: created folder ' folder])
                else
                    disp(['Setup: using existing folder: ' folder])
                end
            end
        end
    end
    clear folder
end

disp('Setup: Path setup done.')
end

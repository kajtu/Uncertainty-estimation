function [minmaxSim] = findSimulationUncertainty(experimentNames,parameters,constants,inds,data,options,bestParams,saveResults,loadResults,expName)
% Load or simulate the prediciton uncertainty from the given parameter
% sets.
%parameters: cell with one or several parameter sets for each patient
%constants: cell with constant values for each patient

    basefolder = split(pwd,'Uncertainty-estimation');
    basefolder = fullfile(basefolder{1},'Uncertainty-estimation');

%% Load saved results if possible
if nargin > 9
    loadfiles=dir(fullfile(basefolder,'Parameters',['minmaxsims' expName '_*']));
else
    loadfiles=dir(fullfile(basefolder,'Parameters','minmaxsims_*'));
end
if loadResults && ~isempty(loadfiles)
    [~,latestfolderInd] = max([loadfiles.datenum]);
    foldernames = {loadfiles.name};
    folderpaths = {loadfiles.folder};
    loadfilename = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    load(loadfilename,'minmaxSim');
    fprintf('findSimulationUncertainty: Loading results from %s \n',loadfilename)    
else
    %% Simulation settings
    step = 0.001;
    
    %% Setup variables
    minmaxSim.States.min = cell(size(parameters));
    minmaxSim.States.max = cell(size(parameters));
    minmaxSim.States.best = cell(size(parameters));
    
    minmaxSim.Observables.min = cell(size(parameters));
    minmaxSim.Observables.max = cell(size(parameters));
    minmaxSim.Observables.best = cell(size(parameters));
    
    minmaxSim.Observables.all = cell(size(parameters));
    minmaxSim.States.all = cell(size(parameters));
    
    minmaxSim.pressurenames = {'pLA','pLV','pAo','pPV'};
    minmaxSim.pressureindexes = [17,18,1,13];
    for pr = 1:length(minmaxSim.pressurenames)
        minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).best = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).best = cell(size(parameters));
    end
    
    %% Run simulations
    w = warning('off','all');
    for p = 1:size(parameters,2)%for each subject
        thisdata = data.(experimentNames{p});
        T = constants(inds{p}.T,p);
        simtime = sort([thisdata.MV.time',0:step:T]);
        simtime = unique(simtime);
        options.x0 = thisdata.IC;
        doSimulation = str2func(thisdata.meta.simulationFunction);
        modelName = thisdata.meta.modelName;

        parameters{p} = parameters{p}';
        %simulate best paramset
        if ~isempty(bestParams)
            simLast = doSimulation(bestParams(:,p),constants(:,p),options,inds{p},simtime,modelName);
            [estimatedSBP,estimatedDBP] = brachialpressure(simLast,inds{p},thisdata);
            minmaxSim.States.best{p} = simLast.x;
            minmaxSim.Observables.best{p} = simLast.y;
            minmaxSim.SBP.best{p} = estimatedSBP;
            minmaxSim.DBP.best{p} = estimatedDBP;
            SVsim=trapz(simLast.t,simLast.x(:,inds{p}.AV));
            EFsim = 100*SVsim/max(simLast.x(:,5));
            minmaxSim.SV.best{p} = SVsim;
            minmaxSim.EF.best{p} = EFsim;
            
            minmaxSim.States.min{p} = simLast.x;
            minmaxSim.States.max{p} = simLast.x;
            minmaxSim.Observables.min{p} = simLast.y;
            minmaxSim.Observables.max{p} = simLast.y;
            minmaxSim.time{p} = simLast.t;
            
            minmaxSim.SBP.max{p} = estimatedSBP;
            minmaxSim.SBP.min{p} = estimatedSBP;
            minmaxSim.DBP.max{p} = estimatedDBP;
            minmaxSim.DBP.min{p} = estimatedDBP;
            minmaxSim.SV.max{p} = SVsim;
            minmaxSim.SV.min{p} = SVsim;
            minmaxSim.EF.max{p} = EFsim;
            minmaxSim.EF.min{p} = EFsim;
            
            for pr = 1:length(minmaxSim.pressurenames)
                thispressure = simLast.y(:,minmaxSim.pressureindexes(pr));
                minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max{p} = mean(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min{p} = mean(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).best{p} = mean(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max{p} = max(thispressure)-min(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min{p} = max(thispressure)-min(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).best{p} = max(thispressure)-min(thispressure);
            end
            
        end
        
        %reduce number of simulated params if they are too many
        numparams=size(parameters{p},2);
        if numparams>2000
            %select paramsets with min & max for each parameter
            [minvals,minind] = min(parameters{p},[],2);
            [maxvals,maxind] = max(parameters{p},[],2);
            %sample the rest eith equal distance to get 2000 params
            allinds = 1:numparams;
            allinds([maxind; minind; numparams]) = []; %also always include last paramset
            indrest=round(linspace(1,length(allinds),2000-length(allinds)));
            indtosim = [maxind; minind; numparams; allinds(indrest)'];
            %always include last paramset
            paramsToSim = parameters{p}(:,indtosim);

            fprintf('The number of simulated parameter sets is reduced from %d to %d\n',numparams,2000)
            % ind=round(linspace(1,numparams,2000));

        else
            paramsToSim = parameters{p};
        end
        
        % Simulate all parameter sets
        fprintf('Simulating %d parameter sets...\n',size(paramsToSim,2))
        for s = 1:size(paramsToSim,2)%for each parameter set
            simLast = doSimulation(paramsToSim(:,s),constants(:,p),options,inds{p},simtime,modelName);
            [estimatedSBP,estimatedDBP] = brachialpressure(simLast,inds{p},thisdata);
            SVsim=trapz(simLast.t,simLast.x(:,inds{p}.AV));
            EFsim = 100*SVsim/max(simLast.x(:,5));
            if isempty(bestParams) && s == 1
                minmaxSim.States.min{p} = simLast.x;
                minmaxSim.States.max{p} = simLast.x;
                minmaxSim.Observables.min{p} = simLast.y;
                minmaxSim.Observables.max{p} = simLast.y;
                minmaxSim.time{p} = simLast.t;
                minmaxSim.SBP.max{p} = estimatedSBP;
                minmaxSim.SBP.min{p} = estimatedSBP;
                minmaxSim.DBP.max{p} = estimatedDBP;
                minmaxSim.DBP.min{p} = estimatedDBP;
                minmaxSim.SV.max{p} = SVsim;
                minmaxSim.SV.min{p} = SVsim;
                minmaxSim.EF.max{p} = EFsim;
                minmaxSim.EF.min{p} = EFsim;
                
                minmaxSim.States.best{p} = NaN(size(simLast.x));
                minmaxSim.Observables.best{p} = NaN(size(simLast.y));
                minmaxSim.SBP.best{p} = NaN;
                minmaxSim.DBP.best{p} = NaN;
                minmaxSim.SV.best{p} = NaN;
                minmaxSim.EF.best{p} = NaN;
                
                for pr = 1:length(minmaxSim.pressurenames)
                    thispressure = simLast.y(:,minmaxSim.pressureindexes(pr));
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max{p} = mean(thispressure);
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min{p} = mean(thispressure);
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).best{p} = NaN;
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max{p} = max(thispressure)-min(thispressure);
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min{p} = max(thispressure)-min(thispressure);
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).best{p} = NaN;
                end
                
                for obs = 1:size(simLast.y,2)
                    minmaxSim.all{obs}{p} =  zeros(size(parameters{p},1),length(simtime));
                end
            else
                % Save the maximum and minimum values for each simulated variable in a struct
                minmaxSim.States.min{p} = min(simLast.x,minmaxSim.States.min{p});
                minmaxSim.States.max{p} = max(simLast.x,minmaxSim.States.max{p});
                minmaxSim.Observables.min{p} = min(simLast.y,minmaxSim.Observables.min{p});
                minmaxSim.Observables.max{p} = max(simLast.y,minmaxSim.Observables.max{p});
                minmaxSim.SBP.max{p} = max(estimatedSBP,minmaxSim.SBP.max{p});
                minmaxSim.SBP.min{p} = min(estimatedSBP,minmaxSim.SBP.min{p});
                minmaxSim.DBP.max{p} = max(estimatedDBP,minmaxSim.DBP.max{p});
                minmaxSim.DBP.min{p} = min(estimatedDBP,minmaxSim.DBP.min{p});
                minmaxSim.SV.max{p} = max(SVsim,minmaxSim.SV.max{p});
                minmaxSim.SV.min{p} = min(SVsim,minmaxSim.SV.min{p});
                minmaxSim.EF.max{p} = max(EFsim,minmaxSim.EF.max{p});
                minmaxSim.EF.min{p} = min(EFsim,minmaxSim.EF.min{p});
                
                for pr = 1:length(minmaxSim.pressurenames)
                    thispressure = simLast.y(:,minmaxSim.pressureindexes(pr));
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max{p} = max(mean(thispressure),minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max{p});
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min{p} = min(mean(thispressure),minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min{p});
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max{p} = max(max(thispressure)-min(thispressure),minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max{p});
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min{p} = min(max(thispressure)-min(thispressure),minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min{p});
                end
            end
        end
    end
    w = warning('on','all');
    if saveResults && nargin > 9
        save(sprintf('%s/Parameters/minmaxsims%s_%s',basefolder,expName,datestr(now,'yymmdd-hhMM')),'minmaxSim','experimentNames');
    elseif saveResults
        save(sprintf('%s/Parameters/minmaxsims_%s',basefolder,datestr(now,'yymmdd-hhMM')),'minmaxSim','experimentNames');
    end
end

end

function [estimatedSBP,estimatedDBP] = brachialpressure(simLast,ind,data)
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
end

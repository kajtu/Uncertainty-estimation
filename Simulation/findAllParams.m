function [expNamesSims,bestparams,bestcosts,meanparams,allokParams,...
    allokCosts,minValuesparams,maxValuesparams,numsincluded,data,inds,constants,medianparams,all,okChi2] = findAllParams(resultsFolder,experimentNames,paramNames,data,inds,constants,expNameBase)
%% Load resulting parameter sets from ESS, PL, and MCMC.

%% Load saved results if possible
%load saved results if it exists
loadResults = 1;
if nargin >6
    loadfiles=dir(fullfile(resultsFolder, ['allokparams' expNameBase '_*']));
else
    loadfiles=dir(fullfile(resultsFolder, 'allokparams_*'));
end
if loadResults && ~isempty(loadfiles)
    % Always load the latest results
    [~,latestfolderInd] = max([loadfiles.datenum]);
    foldernames = {loadfiles.name};
    folderpaths = {loadfiles.folder};
    loadfilename = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    load(loadfilename,'expNamesSims','bestparams','bestcosts','meanparams','medianparams','allokParams',...
        'allokCosts','minValuesparams','maxValuesparams','numsincluded','data','inds','constants','all','okChi2');
    fprintf('findAllParams: Loading existing results from %s \n',loadfilename)
    
    % In patient-specific predictions (eg. when running createFigures_hypertensionT2D_predictions),
    % only the results for one subject at a time is wanted.
    % Thus, we remove all other subjects:
    removenow = ~ismember(expNamesSims,experimentNames);
    if sum(removenow) > 0
        disp('findAllParams: removing all subjects not wanted. Wanted subjects:')
        disp(experimentNames)
        allokParams(removenow) = [];
        allokCosts(removenow) = [];
        bestparams(:,removenow) = [];
        bestcosts(removenow) = [];
        meanparams(:,removenow) = [];
        medianparams(:,removenow) = [];
        constants(:,removenow) = [];
        minValuesparams(:,removenow) = [];
        maxValuesparams(:,removenow) = [];
        inds(removenow) = [];
        numsincluded.MCMC(removenow) = [];
        numsincluded.ESS(removenow) = [];
        all(removenow) = [];
        data = rmfield(data,expNamesSims(removenow));
        expNamesSims(removenow) = [];
    end
        
else
    %pre-define vectors
    expNamesSims = [];
    maxValuesparams = zeros(length(paramNames),length(experimentNames));
    minValuesparams= zeros(length(paramNames),length(experimentNames));
    allokParams = cell(1,length(experimentNames));
    allokCosts = cell(1,length(experimentNames));
    meanparams = zeros(length(paramNames),length(experimentNames));
    medianparams = zeros(length(paramNames),length(experimentNames));
    bestparams= zeros(length(paramNames),length(experimentNames));
    bestcosts= zeros(length(experimentNames),1);
    numsincluded.MCMC = nan(1,length(experimentNames));
    numsincluded.ESS = nan(1,length(experimentNames));
    all = cell(1,length(experimentNames));
    okChi2 = nan(1,length(experimentNames));

    % Load parameters for each experiment e
    for e = 1:length(experimentNames)
        experiment = experimentNames{e};
        disp(experiment)
        %% Load ESS params
        [foundESS,bestessParam,bestessCost,allparamsESS,allcostsESSChi2,bestcostESSchi2] = findESSparams(experiment,resultsFolder,paramNames);

        %% Load MCMC params (not used in this paper)
        foundMCMC=0;
        bestMCMCcost=NaN;

        %% Load PL params
        loadParameters=1;
        [foundPL,allPL,allparamsPL,allcostsPL,bestPLParam,bestPLCost] = findPLparams(experiment,resultsFolder,length(paramNames),loadParameters,paramNames);

        %% best of mcmc, ess, and  PL
        bestparams(:,e) = bestessParam';
        bestcosts(e) = bestcostESSchi2;
        if bestMCMCcost < bestcosts(e)
            bestparams(:,e) = bestMCMCparams;
            bestcosts(e) = bestMCMCcost;
        end
        if bestPLCost < bestcosts(e)
            bestparams(:,e) = bestPLParam';
            bestcosts(e) = bestPLCost;
        end

        %% Combine all params, include the ones below chi2 threshold
        if foundESS
            expNamesSims = [expNamesSims,{experimentNames{e}}];
            % Combine all params
            all{e}.params = allparamsESS(:,2:end);
            all{e}.costs = allcostsESSChi2;
            if foundMCMC
                all{e}.params = [all{e}.params;allMCMCparams'];
                all{e}.costs = [all{e}.costs;allMCMCcosts'];
            else
                disp('findAllParams: no MCMC params found')
            end
            if foundPL
                all{e}.params = [allparamsPL;bestparams(:,e)'];
                all{e}.costs = [allcostsPL;bestcosts(e)];
                disp('findAllParams: only using PL params + bestparam, not ESS and MCMC params.')
            else
                disp('findAllParams: no PL params found')
            end

            % sort all params based on cost
            [all{e}.costs,sortindex] = sort(all{e}.costs);
            all{e}.params = all{e}.params(sortindex,:);

            % Check chi2test
            [dgf,~] = degreesOfFreedom(data.(experimentNames{1}),inds{1});
            thresholdChi2 = chi2inv(0.95,dgf);
            if bestcosts(e) < thresholdChi2
                disp('findAllParams: including params under chi2 (opt + 1 dgf).')
                threshold = bestcosts(e) + chi2inv(0.95,1);
                okChi2(e) = 1;
            else
                %include all params within +10% of best cost
                disp('findAllParams: not ok chi2. including params under 1 dgf of best cost.')
                threshold = bestcosts(e) + chi2inv(0.95,1);
                okChi2(e) = 0;
            end
            okinds = all{e}.costs <= threshold;
            allokParams{e} = all{e}.params(okinds,:);
            allokCosts{e} = all{e}.costs(okinds,:);
            
            meanparams(:,e) = mean(allokParams{e});
            medianparams(:,e) = median(allokParams{e});
            
            minValuesparams(:,e) = min(allokParams{e},[],1);
            maxValuesparams(:,e) = max(allokParams{e},[],1);
            
            if foundMCMC
                numsincluded.MCMC(e) = sum(allMCMCcosts <= bestMCMCcost*1.1);
            else
                numsincluded.MCMC(e) = NaN;
            end
            numsincluded.ESS(e) = sum(allparamsESS(:,1) <= bestMCMCcost*1.1);
        else
            all{e}.params = NaN;
            all{e}.costs = NaN;
        end
    end
    
    %% remove subjects where no parameters were loaded
    if isempty(expNamesSims)
        removeindex = 1:length(experimentNames);
    else
        removeindex = ~ismember(experimentNames,expNamesSims);
    end
    allokParams(removeindex) = [];
    allokCosts(removeindex) = [];
    bestparams(:,removeindex) = [];
    bestcosts(removeindex) = [];
    meanparams(:,removeindex) = [];
    medianparams(:,removeindex) = [];
    constants(:,removeindex) = [];
    minValuesparams(:,removeindex) = [];
    maxValuesparams(:,removeindex) = [];
    inds(removeindex) = [];
    numsincluded.MCMC(removeindex) = [];
    numsincluded.ESS(removeindex) = [];
    all(removeindex) = [];
    [data,keptExpNames] = extractData(expNamesSims,[],data);
    
    %% Save results
    if nargin >6
        savefilename = fullfile(resultsFolder,sprintf('allokparams%s_%s',expNameBase,datestr(now,'yymmdd-HHMM')));
    else
        savefilename = fullfile(resultsFolder,sprintf('allokparams_%s',datestr(now,'yymmdd-HHMM')));
    end
    save(savefilename,'expNamesSims','bestparams','bestcosts','meanparams','medianparams','allokParams',...
        'allokCosts','minValuesparams','maxValuesparams','numsincluded','removeindex','data','inds','constants','all','okChi2');
    
end



end
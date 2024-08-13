function [foundPL,allPL,allparamsPL,allcostsPL,bestPLParam,bestPLCost] = findPLparams(experiment,resultsFolder,numParams,loadParameters,paramNames,paramNamesALL)
% Find the latest folder for PL opt and load the parameters for the experiment.

if loadParameters && nargin < 6
    paramNamesALL = paramNames;%only if specifically have different lists you need to send in all parameter names in the model in order to find the correct pind for the pvalue in extractPLparamfromfolder
end

% set everything to NaN in case it is not found
bestPLParam=NaN;bestPLCost=NaN;allparamsPL_p=NaN;allcostsPL_p=NaN;
allPL = struct;

% find latest patientfolder with PL optimizations
foldersPL = fullfile(resultsFolder,'PL',['E_' experiment]);
if isempty(dir(foldersPL))
    disp(['findPLparams: OBS couldnt load PL parameters for ' experiment, ' (no prev folder)'])
    foundPL = false;
else
    foundPL = false;
    for i = 1:length(paramNames)
        paramfiles = dir(sprintf('%s/allcostsLoaded_%s.mat',foldersPL,paramNames{i}));
        combinedcosts = [];
        for f = 1:length(paramfiles)
            load(fullfile(paramfiles(f).folder,paramfiles(f).name),'allcosts')
            combinedcosts = [combinedcosts;allcosts];
        end

        if ~isempty(combinedcosts)
            [~,indsu] = unique(combinedcosts(:,1));
            combinedcostsUnique=combinedcosts(indsu,:);
            combinedcostsUnique = combinedcostsUnique(~isnan(combinedcostsUnique(:,1)),:);
            for j = 1:length(combinedcostsUnique(:,1))
                pval = combinedcostsUnique(j,1);
                pcosts = combinedcosts(combinedcosts(:,1) == pval,2);
                combinedcostsUnique(j,2) = min(pcosts);
            end

            allPL.(paramNames{i}) = combinedcostsUnique;
            foundPL = true;
        end
    end

    if ~foundPL
        disp(['findPLparams: OBS couldnt load PL parameters for ' experiment, ' (no allcosts files in the folder) ' foldersPL])
    end
end


%%
if foundPL && loadParameters %load the best parameter set for each fixed value
    allparamsPL = [];
    allcostsPL = [];
    bestPLCost = 1e29;
    for param = 1:length(paramNamesALL) %paramNames
        d = dir(fullfile(foldersPL,[paramNamesALL{param} '_*']));
        dfolders = d([d(:).isdir]); % remove files
        pvaluefolders = dfolders(~ismember({dfolders(:).name},{'.','..'})); %remove . and ..
        allparamsPL_p = nan(length(pvaluefolders),numParams);
        allcostsPL_p = nan(length(pvaluefolders),1);
        for i = 1:length(pvaluefolders)
            folderName = fullfile(pvaluefolders(i).folder,pvaluefolders(i).name);
            [pval,allcostsPL_p(i),allparamsPL_p(i,:)] = extractPLparamfromfolder(folderName,numParams,paramNamesALL);

            if allparamsPL_p(i,strcmp(paramNamesALL,'Cpvc')) > 100
                disp('findPLparams: Not including Cpvc > 100, unphysiological value!')
                allparamsPL_p(i,:) = nan(size(allparamsPL_p(i,:)));
                allcostsPL_p(i) = NaN;
            end

            if allparamsPL_p(i,strcmp(paramNamesALL,'k_syst_LV')) > 200
                disp('findPLparams: Not including k_syst_LV > 200, unphysiological value!')
                allparamsPL_p(i,:) = nan(size(allparamsPL_p(i,:)));
                allcostsPL_p(i) = NaN;
            end

            if sum(allparamsPL_p(i,:) <= 0) >0
                disp('findPLparams: removing parameter values <=0, unphysiological value!')
                allparamsPL_p(i,:) = nan(size(allparamsPL_p(i,:)));
                allcostsPL_p(i) = NaN;
            end

            if allcostsPL_p(i) < bestPLCost
                bestPLCost = allcostsPL_p(i);
                bestPLParam = allparamsPL_p(i,:);
            end
        end
        allparamsPL = [allparamsPL;allparamsPL_p];
        allcostsPL = [allcostsPL;allcostsPL_p];
    end
    
end

end


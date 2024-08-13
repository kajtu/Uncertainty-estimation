function [bestcostAll,bestparamAll,bestPLCosts] = collectPLparams(foldersPL,paramNames,paramNamesALL,folderESS)
% Function to load the best cost for each parameter step taken in the
% profile likelihood, by searching for all folders containing the
% optimization for all steps for the specific biomarker and subject. All
% steps together with the best cost for each step is saved in the file
% allcostsLoaded_xx.mat.

% paramNamesALL = {'Cpvc' 'Rpu' 'Rpv' 'Lpv' 'Rtot' 'Ctot' 'ELCo' 'Caa' 'Emax_LA' 'Emax_LV' 'Emin_LA' 'Emin_LV' 'Lao' 'Lav'...
%     'Lmv' 'Ppu' 'Rao' 'Rmv' 'k_diast_LA' 'k_diast_LV' 'k_syst_LA'...
%     'k_syst_LV' 'm1_LA' 'm1_LV' 'm2_LA' 'm2_LV' 'onset_LA' 'onset_LV'};

numParams=length(paramNamesALL);
bestcostAll = 1e29;
bestparamAll = nan(size(paramNamesALL));
bestPLCosts = ones(size(paramNamesALL)).*1e29;
  

%% load optimal parameters found with ESS optimization
if nargin > 3
    folders = dir(folderESS);
    if isempty(folders)
        disp('CollectPLparams: ERROR: no ess found')
    else
        [~,latestfolderInd] = max([folders.datenum]);
        foldernames = {folders.name};
        folderpaths = {folders.folder};
        loadfolderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
        [bestParams,optcost,~,loadedconstants,~,~,~,optcostChi2] = findBestParams(loadfolderName,0,length(paramNamesALL));
        fprintf('CollectPLparams: Loaded ess values with cost %0.2f\n',optcostChi2)
        bestcostAll= optcostChi2;
        bestparamAll = bestParams;
    end
end

%% load the parameters and costs from profile likelihood and save in a collected file
for param = 1:length(paramNames)
    d = dir(fullfile(foldersPL,[paramNames{param} '_*']));
    dfolders = d([d(:).isdir]); % remove files
    pvaluefolders = dfolders(~ismember({dfolders(:).name},{'.','..'})); %remove . and ..
    bestPLCost = 1e29;
    allparamsPL_p = nan(length(pvaluefolders),numParams);

    % add the ESS parameter set
    if bestcostAll ~= 1e29
        allcostsPL_p = nan(length(pvaluefolders)+1,1);
        pvals = nan(length(pvaluefolders)+1,1);
        pind = strcmp(paramNamesALL,paramNames{param});
        allcostsPL_p(end) = bestcostAll;
        pvals(end) = bestparamAll(pind);
    else
        allcostsPL_p = nan(length(pvaluefolders),1);
        pvals = nan(length(pvaluefolders),1);
    end

    % add all PL params
    for i = 1:length(pvaluefolders)
        folderName = fullfile(pvaluefolders(i).folder,pvaluefolders(i).name);
        [pvals(i),allcostsPL_p(i),allparamsPL_p(i,:)] = extractPLparamfromfolder(folderName,numParams,paramNamesALL);

        if strcmp(paramNames{param},'Cpvc') && pvals(i) > 100
            disp('CollectPLparams: Not including Cpvc > 100, unphysiological value!')
            allcostsPL_p(i) = NaN;
            pvals(i) = NaN;
            allparamsPL_p(i,:) = nan(size(allparamsPL_p(i,:)));
        end

        if allparamsPL_p(i,strcmp(paramNamesALL,'k_syst_LV')) > 200
            disp('CollectPLparams: Not including ksystLV > 200, unphysiological value!')
            allcostsPL_p(i) = NaN;
            pvals(i) = NaN;
            allparamsPL_p(i,:) = nan(size(allparamsPL_p(i,:)));
        end

        if pvals(i) <= 0
            disp(['CollectPLparams: removing parameter values <=0 for parameter ' paramNames{param} ', unphysiological value!'])
            allcostsPL_p(i) = NaN;
            pvals(i) = NaN;
            allparamsPL_p(i,:) = nan(size(allparamsPL_p(i,:)));
        elseif allcostsPL_p(i) < bestPLCosts(param)
            bestPLCosts(param) = allcostsPL_p(i);
            bestPLParam = allparamsPL_p(i,:);
        end
    end
    naninds  = isnan(allcostsPL_p);
    allcosts=[pvals(~naninds),allcostsPL_p(~naninds)];
    save(sprintf('%s/allcostsLoaded_%s.mat',foldersPL,paramNames{param}),'allcosts');

    if bestPLCosts(param) < bestcostAll
        bestcostAll = bestPLCosts(param);
        bestparamAll = bestPLParam;
    end
end


%% add the best overall parameter set to all other params (optional)
disp('Adding best overall parameter set to all other biomarker folders (optional)...')
for param = 1:length(paramNames)
    try
        load(sprintf('%s/allcostsLoaded_%s.mat',foldersPL,paramNames{param}),'allcosts');
        allcosts = [allcosts;[bestparamAll(strcmp(paramNamesALL,paramNames{param})),bestcostAll]];
        save(sprintf('%s/allcostsLoaded_%s.mat',foldersPL,paramNames{param}),'allcosts');
    catch
        sprintf('CollectPLparams: Could not find %s/allcostsLoaded_%s.mat to add the best overall paramset to allcosts\n',foldersPL,paramNames{param})
    end
end

end
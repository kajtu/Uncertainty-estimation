function [bestparams,bestcostAll,allparams,bestconstants,meanparams,medianparams,allCostsChi2,bestcostAllchi2] = findBestParams(folderName,loadsummaryfile,numParams,loadMean)
% Load parameters with lowest cost from the given folder called folderName.

if nargin <3
    numParams = 22;
else
    numParams = numParams+1;
end

meanparams = NaN;
medianparams = NaN;
files = dir(fullfile(folderName,'opt-*'));
summaryfile = dir(fullfile(folderName,'bestparam*'));

% some old files needed to be excluded
if contains(folderName,'PL/') || contains(folderName,'PL\')
   files = files([files.datenum] >= 7.391936896527777e+05);
   removedOldFiles = 1;
else
    removedOldFiles = 0;
end

if isempty(files)
    if removedOldFiles
        disp(['findBestParams: ERROR - emtpy load folder (removed old files) ' folderName])
    else
        disp(['findBestParams: ERROR - emtpy load folder ' folderName])
    end
    bestparams  =NaN;
    bestcostAll = NaN;
    allparams  =NaN;
    bestconstants = NaN;
    allCostsChi2=NaN;
    bestcostAllchi2 = NaN;
elseif ~loadsummaryfile || (loadsummaryfile && isempty(summaryfile)) || loadMean
    bestcostAll = 1e100;
    allparams = zeros(length(files),numParams);
    allCostsChi2 = zeros(length(files),1);
    for f = 1:length(files)
        costChi2 = NaN;
        try
            load(fullfile(files(f).folder,files(f).name),'optParam','bestcost','constants','costChi2')
        catch
            optParam = nan(size(allparams(f,2:end)));
            bestcost = 1e99;
            constants = NaN;
            costChi2 = 1e99;
        end
        if isnan(costChi2)
            disp('no chi2cost found')
        end
        if isnan(bestcost)
            disp(fullfile(files(f).folder,files(f).name))
            disp('no bestcost found')
        end
        allparams(f,1) = bestcost;
        allCostsChi2(f) = costChi2;
        if f == 1 && length(optParam) ~= length(allparams(f,2:end)) %if wrong param length
            allparams = zeros(length(files),length(optParam)+1);
        end
        allparams(f,2:end) = optParam;
        if bestcost < bestcostAll
            bestparams = optParam;
            bestcostAll = bestcost;
            bestconstants = constants;
            bestcostAllchi2 = costChi2;
        end
    end
    
    okinds = allparams(:,1) <= bestcostAll*1.10;
    meanparams = mean(allparams(okinds,2:end));
    medianparams = median(allparams(okinds,2:end));
    
    if loadsummaryfile
        save(sprintf('%s/bestparam(%0.3f).mat',files(f).folder,bestcostAll),'bestparams','bestcostAll','bestconstants','allparams','costChi2','bestcostAllchi2')
    end
else
    %load summaryfile
    fprintf('Findbestparams: Loaded summaryfile %s \n',summaryfile(1).name)
    load(fullfile(summaryfile(1).folder,summaryfile(1).name),'bestparams','bestcostAll','bestconstants','allparams')
end
end
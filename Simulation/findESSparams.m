function [foundESS,bestessParam,bestessCost,allparamsESS,allcostsESSChi2,bestcostESSchi2] = findESSparams(experiment,resultsFolder,paramNames)
% Find the latest folder for ESS opt and load the parameters for the experiment.
loadsummaryfile=0;loadmean=0;
% set everything to NaN in case it is not found
bestessParam=NaN;bestessCost=NaN;allparamsESS=NaN;allcostsESSChi2=NaN;bestcostESSchi2=NaN;

% find latest patientfolder with ESS optimizations
loadESSresultsfolder = fullfile(resultsFolder,'ESS',['E_' experiment, '_*']);
foldersESS = dir(loadESSresultsfolder);
if isempty(foldersESS)
    disp(['findAllParams: OBS couldnt load ESS parameters for ' experiment, ' (no prev folder)'])
    foundESS = false;
else
    foundESS = true;
    [~,latestfolderInd] = max([foldersESS.datenum]);
    foldernames = {foldersESS.name};
    folderpaths = {foldersESS.folder};

    folderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    if isempty(dir([folderName,'/opt-*']))
        disp(['findAllParams: OBS couldnt load ESS parameters for ' experiment, ' (no opt files in prev folder)'])
        foundESS = false;
    else
        [bestessParam,bestessCost,allparamsESS,~,~,~,allcostsESSChi2,bestcostESSchi2] = findBestParams(folderName,loadsummaryfile,length(paramNames),loadmean);
        if length(paramNames) == length(bestessParam)+4 %add corr parameters
            bestessParam = [bestessParam 40 40 40 40];
            allparamsESS = [allparamsESS 40+zeros(size(allparamsESS,1),4)];
        end
    end
end

end
function createSampledMeasurementError(numberOfSamples,experimentName)
% experimentName = 'simRandomSystematic' or  simRandomOnly
rng("default") %for reproducibility

basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');
addpath(genpath(basefolder))

savefolder = fullfile(basefolder,'Data',['SampledErrors_' experimentName]);
mkdir(savefolder)

load('dataSimulated.mat','dataSimulated')
trueData = dataSimulated.(experimentName);

for n = 1:numberOfSamples
    disp(100*n/numberOfSamples)
    [~,e] = sampleData(trueData);
    % [~,e] = sampleData(trueData,[],n); %load RR error
    save(fullfile(savefolder,sprintf('sampledMeasurementError%d.mat',n)),'e')
end

fprintf('%d samples created and saved in %s\n',numberOfSamples,savefolder)
end
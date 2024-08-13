function [sampledData,e] = loadSampledMeasurementError(n,experimentName)
%n: which number of the sampled errors with the same setup to load
%experimentName: simRandomSystematic (new method) or simRandomOnly
%(classical method)

% set up paths
basefolder = split(pwd,'Uncertainty-estimation');
basefolder = fullfile(basefolder{1},'Uncertainty-estimation');

%load sampled error
loadfolder = fullfile(basefolder,'Data',['SampledErrors_' experimentName]);
load(fullfile(loadfolder,sprintf('sampledMeasurementError%d.mat',n)),'e')

% load true data
load(fullfile(basefolder,'Data','dataSimulated.mat'),'dataSimulated')
trueData = dataSimulated.(experimentName);

%re-sample the data with the error
[sampledData,e] = sampleData(trueData,e);

end
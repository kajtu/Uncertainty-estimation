function [error] = calcRRsmoothingerror(data)
% Empirical function to estimate the smoothing error due to RR variability
% in 4D flow MRI data. The error is calculated based on the blood flow curves 
% in MV, PV, AV and AA in the given dataset. This is used to estimate data
% uncertainty when not knowing the true error.
% This function is based on 40 data points per flow curve.

datacontent = fieldnames(data);
error = struct;
%% MV
if ismember('MV',datacontent)
    tstep = max(1,round(length(data.MV.time)/40));
    diastdata = data.MV.mean(data.indtdiast:end);
    [~,peakinds] = findpeaks(diastdata,'NPeaks',2,'MinPeakProminence',2,'MinPeakDistance',4);
    peakinds = peakinds+data.indtdiast-1;
    calculatedDiffMV = zeros(size(data.MV.mean));
    calculatedDiffMV(data.indtdiast:end) = -0.10.*data.MV.mean(data.indtdiast:end);

    % The E and A waves are affected the most
    if ~isempty(peakinds) %E-wave is affected more than the A-wave (negatively)
        calculatedDiffMV(peakinds(1)-2*tstep:peakinds(1)+3*tstep) = -0.15.*data.MV.mean(peakinds(1)-2*tstep:peakinds(1)+3*tstep);
    end
    if length(peakinds) > 1 % A-wave (positively)
        if peakinds(2) < 38
            calculatedDiffMV(peakinds(2)-2*tstep:peakinds(2)+3*tstep) = 0.13.*data.MV.mean(peakinds(2)-2*tstep:peakinds(2)+3*tstep);
        else
            disp('A peak too close to timevector end, no extra sd added')
        end
    end

    % The start of diastole (pos)
    calculatedDiffMV(data.indtdiast:data.indtdiast+5*tstep) = 15;

    error.MV.mean = calculatedDiffMV;
    error.MV.time =  data.MV.time;
end

%% PV
%only a small smoothing error, approximated with +-5 ml/s
if ismember('PV',datacontent)
    error.PV.mean = 5.*ones(size(data.PV.mean));
    error.PV.time = data.PV.time;

    %find the PV peaks (for a healthy "textbook" flow)
    [pks,locs] = findpeaks(data.PV.mean,'Npeaks',2);
    Sind = locs(1);
    Dind = locs(2);
    [pks,locs] = findpeaks(-data.PV.mean,'Npeaks',2);
    Aind = locs(2);
    midSD = locs(1);

    tstep = max(1,round(length(data.MV.time)/40));

    %positive error around start and in the end after the A wave (leave as is)
    
    % negative error around S, D, A
    error.PV.mean(Sind-2*tstep:Sind+4*tstep) =  -error.PV.mean(Sind-2*tstep:Sind+4*tstep);
    error.PV.mean(Dind-2*tstep:Aind+2*tstep) =  -error.PV.mean(Dind-2*tstep:Aind+2*tstep);

    % positive error in between S % D (leave as is)
end

%% AV
if ismember('AV',datacontent)
    tstep = max(1,round(length(data.AV.time)/40));
    calculatedDiffAV = zeros(size(data.AV.mean));
    calculatedDiffAV(1:data.indtdiast) = -0.03.*data.AV.mean(1:data.indtdiast);% -3% in general, more errors at peak
    calculatedDiffAV(data.indtdiast:data.indtdiast+2*tstep) = 10;%positive error
    error.AV.mean = calculatedDiffAV;
    error.AV.time =  data.AV.time;
end

%% AA
if ismember('AA',datacontent)
    tstep = max(1,round(length(data.AA.time)/40));
    calculatedDiffAA = zeros(size(data.AA.mean));
    calculatedDiffAA(1:data.indtdiast) = -0.03.*data.AA.mean(1:data.indtdiast);%-3 in general
    calculatedDiffAA(data.indtdiast:data.indtdiast+2*tstep) = 10;

    error.AA.mean = calculatedDiffAA;
    error.AA.time =  data.AA.time;
end

end
function [dgf,degreesoffreedom] = degreesOfFreedom(data,ind)
% Calculate the degrees of freedom in the data (total number of datapoints)

indsyst = 1:ind.tdiast; %syst
inddiast = ind.tdiast:length(data.MV.time); %diast

datanames = {'MV','AV','AC','PV'};
dataparams = {'Emax_LV','Caa','Ctot','Rtot','ELCo'};
degreesoffreedom = ones(1,length(datanames)+2+length(dataparams));%+2 for sbp and dbp

for i = 1:length(datanames)
    if strcmp('MV',datanames{i})
        simind = inddiast;
    elseif strcmp('PV',datanames{i})
        simind = 1:length(data.PV.time);
    else
        simind = indsyst;
    end
    degreesoffreedom(i) = length(simind);
end

dgf = sum(degreesoffreedom);

end
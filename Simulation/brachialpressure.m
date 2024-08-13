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

function sol = simulate_avatar_correction(theta,constants,options,ind,simtimeOrig,modelName)
% Wrapper simulation function used when the model is estimating the systematic errors, instead
% of including them in the estimation data uncertainty.

% add RR smoothing for 5 heartbeats
doPlot=0;
rrNumber = 5; %prev: 10
randomDistribution = 0;
[sol,~,~] = simulate_avatar_RR(theta,constants,options,ind,simtimeOrig,modelName,rrNumber,randomDistribution,doPlot);%simulate_avatar_RR

 % add corrections for other systematic errors
 % offset correction parameters
%-40 since the bounds are -40 to 40 but opt bounds are 0 to 80 = -40+40 to 40+40.
 theta(ind.mvCorr:ind.pvCorr) = theta(ind.mvCorr:ind.pvCorr)-40;

 sol.y(:,ind.mvCorrVol) =  sol.y(:,ind.mvCorrVol) + theta(ind.mvCorr);
 sol.y(:,ind.avCorrVol) =  sol.y(:,ind.avCorrVol) + theta(ind.avCorr);
 sol.y(:,ind.aaCorrVol) =  sol.y(:,ind.aaCorrVol) + theta(ind.aaCorr);
 sol.y(:,ind.pvCorrVol) =  sol.y(:,ind.pvCorrVol) + theta(ind.pvCorr);

end





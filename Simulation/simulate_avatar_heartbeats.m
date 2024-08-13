function [simAllHeartbeats,simLastHeartbeat] = simulate_avatar_heartbeats(theta,constants,options,numHeartBeats,ind,simtime,modelName)
% Do a simulation of the avatar model for a number of heartbeats.
% %First, set parameter values, then
% simulate a number of heartbeats one at a time.
%step = 0.001; %good resolution when plotting simulation
%step = T/39; %comparing to data (since there are 40 timeframes)
%ind is a struct containing indexes for the parameters
%options contains simulation options
T = constants(ind.T);

% create a function handle to the simulation function
if nargin < 7
    modelName = 'avatar_corr';
end
doSimulation = str2func(['simulate_' modelName]);

% Calc normalizing factors for elastance function based on the parameters
norm_factor_LV = calc_norm_factor(T,theta(ind.k_syst_LV),theta(ind.k_diast_LV),theta(ind.m1_LV),theta(ind.m2_LV));
norm_factor_LA= calc_norm_factor(T,theta(ind.k_syst_LA),theta(ind.k_diast_LA),theta(ind.m1_LA),theta(ind.m2_LA));
constants(end) = norm_factor_LV;
constants(end-1) = norm_factor_LA;

%onset_LV: range set to 1-2 instead of -0.5 to 0.5 --> take onset_LV-1.5
theta(ind.onset_LV) = theta(ind.onset_LV) - 1.5;
%onset LA: adapted after onset LV to be close enough
theta(ind.onset_LA) = 1 + theta(ind.onset_LV) - theta(ind.onset_LA);

% offset correction parameters
%-40 since the bounds are -40 to 40 but opt bounds are 0 to 80 = -40+40 to 40+40.
try
    theta(ind.mvCorr:ind.pvCorr) = theta(ind.mvCorr:ind.pvCorr)-40;
catch
    theta =theta;
end

% First simulation of one heartbeat - to get the sizes of x and y
simtime = simtime+T;
options.tstart = T;
simLastHeartbeat = doSimulation(simtime,theta, constants, [], options);%simulate_avatar_HEALTH doSimulation

options.x0 = simLastHeartbeat.x(end,:)';
simlen = length(simLastHeartbeat.t)-1;
simAllHeartbeats.x = zeros(simlen*numHeartBeats,size(simLastHeartbeat.x,2));
simAllHeartbeats.y = zeros(simlen*numHeartBeats,size(simLastHeartbeat.y,2));
simAllHeartbeats.t = zeros(simlen*numHeartBeats,size(simLastHeartbeat.t,2));

simAllHeartbeats.x(1:length(simLastHeartbeat.t),:) = simLastHeartbeat.x;
simAllHeartbeats.y(1:length(simLastHeartbeat.t),:) = simLastHeartbeat.y;
simAllHeartbeats.t(1:length(simLastHeartbeat.t)) = simLastHeartbeat.t;

% The rest of the simulations
for i = 1:numHeartBeats-1
    simLastHeartbeat = doSimulation(simtime,theta, constants, [], options);%simulate_avatar_HEALTH

    %start next sim with end values from this sim
    options.x0 = simLastHeartbeat.x(end,:)';
    
    %save this sim together with the rest
    simAllHeartbeats.x(simlen*i+2:simlen*(i+1)+1,:) = simLastHeartbeat.x(2:end,:);
    simAllHeartbeats.y(simlen*i+2:simlen*(i+1)+1,:) = simLastHeartbeat.y(2:end,:);
    simAllHeartbeats.t(simlen*i+2:simlen*(i+1)+1)   = simLastHeartbeat.t(2:end)+T*i;
end

 % To compare with data we use the last heartbeat when a "steady state" is established.
 % We need to round the time to 3 digits, and start at t=0 to be able to
 % compare with data.
simLastHeartbeat.t = round(simLastHeartbeat.t-T,3);
simAllHeartbeats.t = simAllHeartbeats.t - T;
end
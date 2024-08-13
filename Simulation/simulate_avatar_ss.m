function [sol] = simulate_avatar_ss(theta,constants,options,ind,simtime,modelName)
% Do a simulation of the avatar model until steady state. 
% First, set parameter values, 
% then simulate a number of heartbeats one at a time, until steady state is
% reached. Save the last heartbeat simulation.
%step = 0.001; %good resolution when plotting simulation
%step = T/39; %comparing to data (since there are 40 timeframes)
%ind is a struct containing indexes for the parameters
%options contains simulation options

numHeartBeats = 20; % simulates until steady state or max 20 heartbeats

T = constants(ind.T);

% create a function handle to the simulation function
if nargin < 6
    modelName = 'avatar_corr';
end
doSimulation = str2func(['simulate_' modelName]);

% Calc normalizing factors for elastance function based on the parameters
constants(end) = calc_norm_factor(T,theta(ind.k_syst_LV),theta(ind.k_diast_LV),theta(ind.m1_LV),theta(ind.m2_LV));
constants(end-1) = calc_norm_factor(T,theta(ind.k_syst_LA),theta(ind.k_diast_LA),theta(ind.m1_LA),theta(ind.m2_LA));

%onset_LV: range set to 1-2 instead of -0.5 to 0.5 --> take onset_LV-1.5
theta(ind.onset_LV) = theta(ind.onset_LV) - 1.5;
%onset LA: adapted after onset LV to be close enough
theta(ind.onset_LA) = 1 + theta(ind.onset_LV) - theta(ind.onset_LA);

% offset correction parameters
%-40 since the bounds are -40 to 40 but opt bounds are 0 to 80 = -40+40 to 40+40.
try
    theta(ind.mvCorr:ind.pvCorr) = theta(ind.mvCorr:ind.pvCorr)-40;
catch
    theta = theta;
end

simtime = simtime+T;
options.tstart = T;

sol.x = NaN;
for i = 1:numHeartBeats
    solbefore = sol.x;
    sol = doSimulation(simtime,theta, constants, [], options);
    
    %if reaches steady state earlier, there is no need to do extra
    %simulations
    if abs(sum(solbefore(:) - sol.x(:))) < 1 %abs(sum(sol.x(end,:) - options.x0')) < 0.001
        break
    end
    %start next sim with end values from this sim
    options.x0 = sol.x(end,:)';
end

 % To compare with data we use the last heartbeat when a "steady state" is established.
 % We need to start at t=0 to be able to compare with data.
sol.t = sol.t-T;

end
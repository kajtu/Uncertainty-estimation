function [Ctot,Rtot,Emax_LV] = calculateDatabasedParameters(data)
% Calculate parameters on a subject-specific basis based on the input data

%% Emax LV
rho_blood=1.06; %g/ml (constant)
B_av=rho_blood/(2*data.parameters.ELCo.mean^2);
maxQav = max(data.AV.mean);
deltaPconst = 0.06/133.322;
deltaP_av = B_av*(maxQav^2)*deltaPconst; %Pressure gradient accross the aortic valve
k=0.9; %Constant to relate systolic blood pressure and end systolic blood pressure: Pes=k*Ps 
V0_LV=10; %(mL) Mynard et al. (constant)
Emax_LV=k*(data.SBP.mean+deltaP_av)/(data.extra.ESV.mean-V0_LV); %Maximal elastance of the LV


%% Ctot
PP=data.SBP.mean-data.DBP.mean;    
meanQav = mean(data.AV.mean);
Ctot=meanQav/PP;  %Total compliance of the system


%% Rtot
HR=60*1/data.extra.T.mean;  %Heart rate (beats/min)
MAP=data.DBP.mean+(1/3+HR*0.0012)*(data.SBP.mean-data.DBP.mean); % Beléns calculation of MAP, takes HR in consideration azminia M, Trivedi A, Molnar J, Elbzour M, Guerrero M, Salem Y, Ahmed A, Khosla S, Lubell DL . Validation of a new formula for mean arterial pressure calculation: the new formula is superior to the standard formula. Catheter Cardiovasc Interv 2004; 63: 419–425.
Rtot=MAP/meanQav;     %Total resistance of the system


end
% Compile model files into mexfiles
% This automatically creates the files: ami_avatar_corr.mex.. and simulate_avatar_corr.m

clear 
clear mex
 
startup

path = pwd;

amiwrap('avatar_corr','avatar_syms_corr',path); 

disp('Models generated')


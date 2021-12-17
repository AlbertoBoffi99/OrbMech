% ASSIGNMENT 1
%   minimize the delta velocity required to depart from Jupiter and arrive
%   on Venus through a fly-by around the Earth

% TO DO
% - add check on Earth atmosphere
% - add check on rp_norm to be higher than RE

%-------------------------------------------------------------------------%

% initialize MATLAB space
clear all
close all
clc

% clock started
tic
fprintf('Execution started ...\n');

%-------------------------------------------------------------------------%

%% INPUT

% run configuration file, containing all numerical data
run config.m

%% OPTIMIZATION

% declare flyby modelling function as handle
JEV = @(dates) JEV_lamb_fb_lamb(dates,astro,options);
% fmincon guess for for cicle
fmincon_guess = optim.fmincon_guess;

for i = 1:1:optim.n_opt
    
    % optimization of Dv using fmincon
    [results.fmincon_dates(i,:), results.fmincon_Dv(i)] = fmincon(JEV, fmincon_guess, [], [], [], [], [], [], [], options.fmincon_options);
    
    % B constraint matrix for GA 
    optim.ga_Bcon = [ -(results.fmincon_dates(i,1) - optim.ga_SW);
                        (results.fmincon_dates(i,1) + optim.ga_SW);
                        -(results.fmincon_dates(i,2) - optim.ga_SW);
                        (results.fmincon_dates(i,2) + optim.ga_SW);
                        -(results.fmincon_dates(i,3) - optim.ga_SW);
                        (results.fmincon_dates(i,3) + optim.ga_SW)  ];
    
    % optimization of Dv using ga
    [results.ga_dates(i,:), results.ga_Dv(i,:)] = ...
        ga(JEV, optim.ga_nvars, optim.ga_Acon, optim.ga_Bcon, ...
        [], [], [], [], [], options.ga_options);
    % update fmincon guess according to ga result
    fmincon_guess = results.ga_dates(i,:);

end

clear fmincon_guess i

%-------------------------------------------------------------------------%


% clock stopped
fprintf('\nExecution concluded\n')
toc



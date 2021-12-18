% ASSIGNMENT 1
%   minimize the delta velocity required to depart from Jupiter and arrive
%   on Venus through a fly-by around the Earth

% TO DO
% - add check on Earth atmosphere
% - add check on rp_norm to be higher than RE
% - add plot

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
JEV = @(dates) JEV_lamb_fb_lamb(dates, astro, options);

for i = 1:1:10
    i
% optimization of Dv using ga
[results.ga_dates, results.ga_Dv] = ...
    ga(JEV, optim.ga_nvars, optim.ga_Acon, optim.ga_Bcon, ...
    [], [], [], [], [], options.ga_options);

%fprintf('\n Genetic algorithm optimization conlcuded\n')

% update fmincon guess according to ga result
fmincon_guess = results.ga_dates;

% optimization of Dv using fmincon
[results.fmincon_dates, results.fmincon_Dv] = fmincon(JEV, fmincon_guess, [], [], [], [], [], [], [], options.fmincon_options);

results.fmincon_Dv
results.ga_Dv

%fprintf('\n Gradient method optimization conlcuded\n')
end

%-------------------------------------------------------------------------%


% clock stopped
fprintf('\nExecution concluded\n')
toc

%% PLOTTING


if options.plot
    fprintf('\nPlotting started\n')
    run plotting.m
    fprintf('\Plotting concluded\n')
end



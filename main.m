% ASSIGNMENT 1
%   minimize the delta velocity required to depart from Jupiter and arrive
%   on Venus through a fly-by around the Earth

% TO DO
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

for index = 1:1:options.noptim

fprintf('\nIteration number %d started\n', index);

% optimization of Dv using ga
[results.ga_dates(index, :), results.ga_Dv(index)] = ...
    ga(JEV, optim.ga_nvars, optim.ga_Acon, optim.ga_Bcon, ...
    [], [], [], [], [], options.ga_options);

% calling global function results
global results

fprintf('\nGenetic algorithm optimization conlcuded\n')

% update fmincon guess according to ga result
fmincon_guess = results.ga_dates(index,:);

% optimization of Dv using fmincon
[results.fmincon_dates(index, :), results.fmincon_Dv(index)] = fmincon(JEV, fmincon_guess, [], [], [], [], [], [], [], options.fmincon_options);

fprintf('\nGradient method optimization conlcuded\n')
end

clear index

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



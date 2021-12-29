% ASSIGNMENT 2
%   Study repeating and perturbed orbits around the Earth. Compare the
%   results with a real satellite trajectory

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

%% NOMINAL ORBIT

% cartesian initial conditions
[nominal.r0 nominal.v0] = kep2car(nominal.a0, nominal.e0, deg2rad(nominal.i0), deg2rad(nominal.Ome0), deg2rad(nominal.ome0), deg2rad(nominal.f0), astro.muE);

% one period interval of integration 
nominal.tspan = linspace(real.mjd2000_start*24*3600, real.mjd2000_start*24*3600 + nominal.T, 1000);

% integration
[nominal.time nominal.car] = ode113(@(t,s) tbp_ode(t,s,astro.muE), nominal.tspan, [nominal.r0; nominal.v0], options.ode);


%% GROUNDTRACKS CALCULATION

% semi-major axis of an unperturbed repeating groundtrack
repun.a = repeatingGroundTracks(astro.k, astro.m, astro.muE, astro.wE);
% keplerian elements of an unperturbed repeating groundtrack
repun.kep = [repun.a nominal.kep(2:end)];

switch choice.GTperiod
    case 1
        temp.tspan = linspace(real.mjd2000_start*24*3600, (real.mjd2000_start+1)*3600*24, 100000);
    case 2
        temp.tspan = linspace(real.mjd2000_start*24*3600, (real.mjd2000_start+10)*3600*24,100000);
    case 3
        temp.tspan = linspace(real.mjd2000_start*24*3600, real.mjd2000_start*24*3600 + nominal.T, 100000);
end

[nominal.alpha, nominal.delta, nominal.lon, nominal.lat] = groundTrack...
    (nominal.kep, deg2rad(nominal.lon0), temp.tspan, astro.wE, astro.muE, nominal.t0);

[repun.alpha, repun.delta, repun.lon, repun.lat] = groundTrack...
    (repun.kep, deg2rad(nominal.lon0), temp.tspan, astro.wE, astro.muE, nominal.t0);

[per.alpha, per.delta, per.lon, per.lat] = groundTrack_perturbed...
    (nominal.kep, deg2rad(nominal.lon0), temp.tspan, astro.wE, astro.muE, astro.RE, astro.J2, astro.cr, astro.A2m, nominal.t0);

[repper.alpha, repper.delta, repper.lon, repper.lat] = groundTrack_perturbed...
    (repun.kep, deg2rad(nominal.lon0), temp.tspan, astro.wE, astro.muE, astro.RE, astro.J2, astro.cr, astro.A2m, nominal.t0);

clear temp

%% INTEGRATION OF PERTURBED ORBITS (SRP and J2)

% time interval of integration 
temp.tspan_3year = linspace(real.mjd2000_start*24*3600, real.mjd2000_start*24*3600 + 3*365*3600*24, 40000);

% time interval of integration comparison with real satellite
temp.tspan_real = real.day*24*3600; 

% propagation of Cartesian elements
[per.time_car, per.car] = ode113( @(t,s) ode_perturbed ...
    (t, s, astro.muE,  astro.RE, astro.J2, astro.cr, astro.A2m), temp.tspan_3year, [nominal.r0; nominal.v0], options.ode ); 

% propagation of Keplerian elements
[per.time_kep, per.kep_gauss] = ode113( @(t,s) gauss_propagation...
    (t, s, astro.muE, astro.RE, astro.J2, astro.cr, astro.A2m), temp.tspan_3year, nominal.kep', options.ode); 

% propagation of Keplerian elements comparison with real satellite
[per.time_kep_real, per.kep_gauss_real] = ode113(@(t,s) gauss_propagation...
    (t, s, astro.muE, astro.RE, astro.J2, astro.cr, astro.A2m), temp.tspan_real, real.kep0', options.ode); 

% Keplerian elements to cartesian elements
for j = 1:length(temp.tspan_3year)
    [per.car_gauss(j,1:3), per.car_gauss(j, 4:6)] = kep2car(per.kep_gauss(j,1), per.kep_gauss(j,2),...
        per.kep_gauss(j,3), per.kep_gauss(j,4), per.kep_gauss(j,5), per.kep_gauss(j,6), astro.muE);
end

% Cartesian elements to keplerian elements
for k = 1:length(temp.tspan_3year)
    [ per.kep(k,1),  per.kep(k,2),  per.kep(k,3),  per.kep(k,4),  per.kep(k,5),  per.kep(k,6) ] = ...
        car2kep(per.car(k,1:3), per.car(k,4:6), astro.muE);
end


% to unwrap RAAN, argument of perigee and true anomaly for plotting
per.kep_unwr(:,1:3) = per.kep(:,1:3);
per.kep_unwr(:,6) = unwrap(per.kep(:,6));
per.kep_unwr(:,4) = unwrap(per.kep(:,4));
per.kep_unwr(:,5) = unwrap(per.kep(:,5));

%-------------------------------------------------------------------------%

% clock stopped
fprintf('\nExecution concluded\n')
toc

%% PLOTTING

% if the user asks for plotting 
if options.plot

    fprintf('\nPlotting started ...\n')

    % plotting function
    run plotting.m
    
    fprintf('Plotting concluded\n')

end


%% SAVE RESULTS

% if the user asks for saving results 
if options.save

    save('.\Results\workspace.mat');
    fprintf('Results saved\n');

end

% clear workspace
clearvars -except astro choice fig filt nominal options per real repper repun

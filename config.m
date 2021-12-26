% CONFIGURATION FILE
% In this file all numerical data is declared

%-------------------------------------------------------------------------%

%% STRUCUTRES DECLARATION

global results out
% intila minimum delta velocity
results.Dv_min = astroConstants(5);
% initial dates (just not to have NaN as first date)
results.dates = [1e6 1e6 1e6];

%% ASTRONOMICAL DATA

% Earth gravitational constant
astro.muE = astroConstants(13);
% Sun gravitational constant
astro.muS = astroConstants(4);
% Jupiter gravitational constant
astro.muJ = astroConstants(15);
% Venus gravitational constant
astro.muV = astroConstants(12);
% Earth radius
astro.RE = astroConstants(23);
% height of Earth atmosphere
astro.h_atm = 100;

%% MINIMIZATION SETTINGS

% A constraint matrix for GA and fmincon
optim.Acon = [ -1 0 0;
               0 -1 0;
               0 0 -1;
               1 0 0;
               0 1 0;
               0 0 1];

% number of variables for GA
optim.ga_nvars = 3;

% mesh settings:
% number of elements in departure window vector
optim.nmesh_dep = 10;
% number of elements in each flyby window vector
optim.nmesh_fly = 5;
% number of elements in each arrival window vector
optim.nmesh_arr = 10;

% number of GA and fmincon optmization in meshpgafmincon.m
optim.noptim = 10;

%% DATES

% departure window
departure.date_min = date2mjd2000([2025 08 01 00 00 00]); 
departure.date_max = date2mjd2000([2040 12 01 00 00 00]);
% departure semi-window to open GA and fmincon window around a first guess
departure.SW = 0.5*(departure.date_max - departure.date_min)/optim.nmesh_dep;
% departure minimum ToF retreived from 1st arc's pork-chop plot
departure.tof_min = 100;
% departure maximum ToF retreived from 1st arc's pork-chop plot
departure.tof_max = 2500;

% flyby window
fly.date_min = date2mjd2000([ 2030 08 01 00 00 00]);
fly.date_max = date2mjd2000([ 2050 12 01 00 00 00]);
% flyby semi-window to open GA and fmincon window around a first guess
fly.SW = 0.5*(fly.date_max - fly.date_min)/optim.nmesh_fly;

% arrival date
arrival.date_min = date2mjd2000([2040 01 01 00 00 00]);  
arrival.date_max = date2mjd2000([2065 08 01 00 00 00]);
% arrival semi-window to open GA and fmincon window around a first guess
arrival.SW = 0.5*(arrival.date_max - arrival.date_min)/optim.nmesh_arr;
% arrival minimum ToF retreived from 2nd arc's pork-chop plot
arrival.tof_min = 100;
% arrival maximum ToF retreived from 2nd arc's pork-chop plot
arrival.tof_max = 2500;

%% FUNCTION OPTIONS

% random numbre generator settings
rng ('shuffle')
% options for GA
options.ga_options = optimoptions("ga","ConstraintTolerance", 1e-6,"PopulationSize", 150, "MaxGenerations", 30, "FunctionTolerance", 1e-6, "Display", "off"); %'PlotFcn' {@gaplotbestf,@gaplotstopping});
% options for fmincon
options.fmincon_options = optimoptions("fmincon", "Display", "off");
% options for fsolve
options.fsolve_options = optimset("Display", "off");
% options for ode
options.ode_options = odeset ('RelTol', 1e-3, 'AbsTol', 1e-3);
% options for lambertMR
options.orbitType = 0;
options.Nrev = 0;
options.Ncase = 0;
options.LambOptions = 1;
% plotting options
options.plot = 1;

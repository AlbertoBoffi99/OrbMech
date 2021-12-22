% CONFIGURATION FILE
% In this file all numerical data is declared

%-------------------------------------------------------------------------%

%% STRUCUTRES DECLARATION

global results optim
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

% A constraint matrix for GA
optim.ga_Acon = [ -1 0 0;
                  0 -1 0;
                  0 0 -1;
                  1 0 0;
                  0 1 0;
                  0 0 1];

% number of variables for GA
optim.ga_nvars = 3;

% lower bound for GA
optim.ga_lb = [];

% perigee height flag in GA
optim.ga_nanflag = 0;

% mesh settings
optim.nmesh_dep = 100;
optim.nmesh_fly = 25;
optim.nmesh_arr = 100;

% optimizations number
optim.noptim = 10;

%% DATES

% departure window
departure.date_min = date2mjd2000([2025 08 01 00 00 00]); 
departure.date_max = date2mjd2000([2040 12 01 00 00 00]);
departure.SW = 0.5*(departure.date_max - departure.date_min)/optim.nmesh_dep;
departure.tof_min = 100;
departure.tof_max = 2500;

% flyby window
fly.date_min = date2mjd2000([ 2030 08 01 00 00 00]);
fly.date_max = date2mjd2000([ 2050 12 01 00 00 00]);
fly.SW = 0.5*(fly.date_max - fly.date_min)/optim.nmesh_fly;

% arrival date
arrival.date_min = date2mjd2000([2040 01 01 00 00 00]);  
arrival.date_max = date2mjd2000([2065 08 01 00 00 00]);
arrival.SW = 0.5*(arrival.date_max - arrival.date_min)/optim.nmesh_arr;
arrival.tof_min = 100;
arrival.tof_max = 2500;

%% FUNCTION OPTIONS

% random numbre generator settings
% rng default
% options for GA
options.ga_options = optimoptions("ga","ConstraintTolerance", 1e-6,"PopulationSize", 100, "MaxGenerations", 30, "FunctionTolerance", 1e-6, "Display", "off"); % 'PlotFcn', {@gaplotbestf,@gaplotstopping}
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

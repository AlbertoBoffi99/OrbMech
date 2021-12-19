% CONFIGURATION FILE
% In this file all numerical data is declared

%-------------------------------------------------------------------------%

%% STRUCUTRES DECLARATION

global results

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

%% DATES

% departure window
departure.date_min = date2mjd2000([2025 08 01 00 00 00]); 
departure.date_max = date2mjd2000([2040 12 01 00 00 00]); 

% flyby window
fly.date_min = date2mjd2000([ 2030 08 01 00 00 00]);
fly.date_max = date2mjd2000([ 2050 12 01 00 00 00]);

% arrival date
arrival.date_min = date2mjd2000([2040 01 01 00 00 00]);  
arrival.date_max = date2mjd2000([2065 08 01 00 00 00]);

%% MINIMIZATION SETTINGS

% A constraint matrix for GA
optim.ga_Acon = [-1 1 0;
                  1 -1 0;
                  0 -1 1;
                  0 1 -1;
                  -1 0 0;
                  0 -1 0;
                  0 0 -1;
                  1 0 0;
                  0 1 0;
                  0 0 1];

% B constraint matrix for GA 
optim.ga_Bcon = [2500; 0.01; 2000; 0.01; -departure.date_min; -fly.date_min; ...
    -arrival.date_min; departure.date_max; fly.date_max; arrival.date_max];

% number of variables for GA
optim.ga_nvars = 3;


%% FUNCTION OPTIONS

% random numbre generator settings
% rng default
% options for GA
options.ga_options = optimoptions("ga","ConstraintTolerance", 1e-6,"PopulationSize", 100, "MaxGenerations", 30, "FunctionTolerance", 1e-6, "Display", "off", 'PlotFcn',{@gaplotbestf,@gaplotstopping});
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
% optimizations number
options.noptim = 10;
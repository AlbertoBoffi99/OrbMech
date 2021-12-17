% CONFIGURATION FILE
% In this file all numerical data is declared

%-------------------------------------------------------------------------%

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
optim.ga_Acon = [-1 0 0;
                  1 0 0;
                  0 -1 0;
                  0 1 0;
                  0 0 -1;
                  0 0 1;];
% number of variables for GA
optim.ga_nvars = 3;
% semi width of the optimization window for GA around fmincon result
optim.ga_SW = 1000;
% fmincon initial guess
optim.fmincon_guess = [1.3120e4 1.4056e4 1.4610e4];
% number of optimization iterations
optim.n_opt = 5;

%% FUNCTION OPTIONS

% options for GA
options.ga_options = optimoptions("ga","ConstraintTolerance", 1e-6,"PopulationSize", 50, "MaxGenerations", 100, "FunctionTolerance", 1e-6, "Display", "off");
% options for fmincon
options.fmincon_options = optimoptions("fmincon", "Display", "off");
% options for fsolve
options.fsolve_options = optimset("Display", "off");
% options for lambertMR
options.orbitType = 0;
options.Nrev = 0;
options.Ncase = 0;
options.LambOptions = 1;
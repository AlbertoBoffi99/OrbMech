% CONFIGURATION FILE
% In this file all numerical data is declared

%-------------------------------------------------------------------------%

%% ASTRO CONSTANTS and PARAMETERS

% Earth's rotational angular velocity
astro.wE = deg2rad(15.04)/3600; % [rad/s]

% Earth's gravitational parameter
astro.muE = astroConstants(13);

% Earth's radius
astro.RE = astroConstants(23);

% Gravitatonal field constant of the Earth
astro.J2 = astroConstants(9);

% parameters for repeating groundtracks
astro.k = 2;
astro.m = 1;

% SRP parameters
astro.cr = 1; 
astro.A2m = 3.000; % [m^2/kg]

%% ORBIT DATA

% selected date for starting the simulation
real.date_start = [2018 12 20 0 0 0];
% convert it to mjd2000
real.mjd2000_start = date2mjd2000(real.date_start);

% real satellite ephemerides import and save [ATLAS 5 CENTAUR R/B NORAD
% 28473]
real.data = importdata("ephemerides.txt", '\t');

z = 1;
for j = 2 : 2 : length(real.data)
    real.kep(z,3) = str2double(extractBetween(real.data(j),9,16));
    real.kep(z,4) = str2double(extractBetween(real.data(j),18,25));
    real.kep(z,2) = str2double(extractBetween(real.data(j),27,33))*1e-7;
    real.kep(z,5) = str2double(extractBetween(real.data(j),35,42));
    real.kep(z,1) = (astro.muE/(str2double(extractBetween(real.data(j),53,63))*2*pi/(24*3600))^2)^(1/3);
    z = z+1;
end
clear z j

for z = 3 : 5
    real.kep(:,z) = unwrap(real.kep(:,z));
end

clear z

% initial condition for integration propagation vs real satellite
real.kep0 = [real.kep(1,1) real.kep(1,2) deg2rad(real.kep(1,3)) deg2rad(real.kep(1,4)) deg2rad(real.kep(1,5)) 0]; 

% extract time vector from ephemerides
z = 1;
real.sidereal_year = 365 + 6/24 + 9/(60*24) + 9.54/(60*24*60); % [days]
for j = 1 : 2 : length(real.data)
    real.time(z) = str2double(extractBetween(real.data(j),19,32));
    real.year(z) = 2000 + floor(real.time(z)/1000);
    real.day(z) = real.time(z) - (real.year(z)-2000)*1000 + real.sidereal_year*floor(real.time(z)/1000); % [days]
    if z > 1
        if real.day(z) - real.day(z-1) <= 0
            real.day(z) = real.day(z) + 1e-10;
        end
    end
    z = z+1;
end
clear z

% assigned initial condition
nominal.a0 = 2.5952e+4; % [km]
nominal.e0 = 0.5555;
nominal.i0 = 17.7395; % [deg]

% satellite period
nominal.T = 2*pi/sqrt(astro.muE/nominal.a0^3); % [s]

% other Keplerian elements to be supposed 
nominal.Ome0 = real.kep(1,4); % [deg]
nominal.ome0 = real.kep(1,5); % [deg]
nominal.f0 = 30; % [deg]

% collect them into a vector
nominal.kep = [nominal.a0, nominal.e0, deg2rad(nominal.i0), deg2rad(nominal.Ome0), deg2rad(nominal.ome0), deg2rad(nominal.f0)];

% initial longitude 
nominal.lon0 = 0; % [deg]

% initial time 
nominal.t0 = 0; % [s]

% perigee height
nominal.hp = 5164.654; % [km]

%% CHOICES

% Groundtracks integration time interval:
%       1 for 1 day
%       2 for 10 days
%       3 for 1 period
choice.GTperiod = 2;

% Keplerian element to filter:
%       1 for semi-major axis 'a'
%       2 for eccentricity 'e'
choice.filtering1 = 1;

% Keplerian element to filter:
%       1 for inclination 'i'
%       2 for RAAN 'Omega'
%       3 for pericenter anomaly 'omega'
%       4 for true anomaly 'f'
choice.filtering2 = 1;

%% OTHER OPTIONS

% ode solver options
options.ode = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);

% saving options
options.save = 0;

% plotting options
options.plot = 1;

% this script changes all interpreters from tex to latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
clear list_factory index_interpreter default_name
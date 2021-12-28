% selected date for starting the simulation
date_start = [2018 12 20 0 0 0]; 
MJD2000 = date2mjd2000(date_start); % [days]

% Earth's gravitational parameter
muE = astroConstants(13);

% real satellite
A = importdata("RealKeplerian1.txt", '\t');

Kep_real = zeros(length(A)/2, 5);
z = 1;
for j=2:2:length(A)
    time_real(z) = str2double()
    Kep_real(z,3) = str2double(extractBetween(A(j),9,16));
    Kep_real(z,4) = str2double(extractBetween(A(j),18,25));
    Kep_real(z,2) = str2double(extractBetween(A(j),27,33))*1e-7;
    Kep_real(z,5) = str2double(extractBetween(A(j),35,42));
    Kep_real(z,1) = (muE/(str2double(extractBetween(A(j),53,63))*2*pi/(24*3600))^2)^(1/3);
    z = z+1;
end
% 
% z = 1;
% for j=1:2:length(A)
%     time_real(z) = str2double(extractBetween(A(j),19,32));
%     year_real(z) = 2000 + floor(time_real(z)/1000);
%     day_real(z) = time_real(z) - (year_real(z)-2000)*1000;
%     z = z+1;
% end




% Earth's rotational angular velocity
omeE = deg2rad(15.04)/3600; % [rad/s]

% Earth's radius
Re = astroConstants(23);

% Gravitatonal field constant of the Earth
J2 = astroConstants(9);

% parameters for repeating groundtracks
k = 2;
m = 1;

% SRP parameters
cr = 1; 
A2m = 3.000; % [m^2/kg]

%% nominal orbit

% assigned initial condition
a0 = 2.5952e+4; % [km]
e0 = 0.5555;
i0 = 17.7395; % [deg]

% other Keplerian elements to be supposed 
Ome0 = Kep_real(1,4); % [deg]
ome0 = Kep_real(1,5); % [deg]
f0 = 0; % [deg]

Kep_elements = [a0, e0, deg2rad(i0), deg2rad(Ome0), deg2rad(ome0), deg2rad(f0)];

% initial longitude 
lon0 = 0; % [deg]

% initial time 
t0 = 0; % [s]

% perigee height
hp = 5164.654; % [km]

% ode solver options
options_sol = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);
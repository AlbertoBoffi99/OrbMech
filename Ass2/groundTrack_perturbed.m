function [alpha, delta, lon, lat] = groundTrack_perturbed (K_elements, lon0, tspan, omeP, muP, Rp, J2, cr, A2m, t0)
% Propagation of satellite's right ascension, declination, longitude and 
% latitude during its motion considering J2 and SRP perturbations
% 
% PROTOTYPE
%     [alpha, delta, lon, lat] = groundTrack_perturbed (K_elements, lon0, tspan, omeP, muP, Rp, J2, cr, A2m, t0)
%     
% INPUT 
%     K_elements [array]      orbital initial Keplerian elements for the
%                             integration of the satellite orbit
%     lon0 [scalar]         initial satellite longitude [rad]                                 
%     tspan [array]           time instants of orbit integration [s]
%     omeP [scalar]         planet's rotational angular velocity [rad/s]
%     muP [scalar]          gravitational constant of the planet [km^3/s^2]
%     Rp [scalar]           mean radius of the planet [km]
%     J2 [scalar]           gravitatonal field constant of the Earth
%     cr [scalar]           factor dependent upon the optical properties of
%                             the spacecraft surface
%     A2m[scalar]           surface to mass rate [m^2/kg]
%     t0 [scalar]           initial time for ground track calculation [s]
%
% OUTPUT
%     alpha [scalar]           satellite's right ascension [rad]
%     delta [scalar]           satellite's declination [rad]
%     lon [scalar]             satellite's longitude [rad]
%     lat [scalar]             satellite's latitude [rad]
%
% CONTRIBUTORS
%     Alberto Boffi, Enrico Raviola, Andrea Campagna, Luca Ciavirella
% 
% VERSION
%     29-12-2021: v01.0
%-------------------------------------------------------------------------%
%% INITIAL SETTINGS
% from keplerian to cartesian elements
[r0,v0] = kep2car(K_elements(1), K_elements(2), K_elements(3), K_elements(4), K_elements(5), K_elements(6), muP);

C_elements = [r0; v0];

%% ORBITAL ELEMENTS INTEGRATION
% setting ode options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% integration
[ ~, state ] = ode113( @(t,s) ode_perturbed(t,s,muP,Rp,J2,cr,A2m), tspan, C_elements, options );

%% OUTPUTS
r = sqrt((state(:,1)).^2 + (state(:,2)).^2 + (state(:,3)).^2); % [km]

% satellite's declination
delta = asin(state(:,3)./r);

% satellite's right ascension
alpha = [];
for k = 1:length(tspan)
    
    if state(k,2)/r(k) > 0
        alpha(k) = acos((state(k,1)/r(k))/cos(delta(k)));
    else
        alpha(k) = 2*pi - acos((state(k,1)/r(k))/cos(delta(k)));
    end

end 

% satellite's longitude
lon_E = lon0 + omeP*(tspan-t0);

lon = alpha - lon_E;

% satellite's latitude
lat = delta;

end
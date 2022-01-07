function [ds] = ode_perturbed (t, s, muP, Re, J2, cr, A2m)
% Propagation of Cartesian elements considering the effect of J2 and SRP 
% 
% Function to be integrated by the ODE solver in order to propagate 
% Cartesian elements in time
% 
% PROTOTYPE
%     [ds] = ode_perturbed (t, s, muP, Re, J2, cr, A2m)
%     
% INPUT
%     muP [scalar]   gravitational constant of the planet [km^3/s^2]
%     Re [scalar]    mean radius of the planet [km]
%     J2 [scalar]    gravitatonal field constant of the Earth
%     cr [scalar]    factor dependent upon the optical properties of
%                      the spacecraft surface
%     A2m[scalar]    surface to mass rate [m^2/kg]
%
% OUTPUT
%     ds [6x1]           Cartesian elements differential equations
%
% CONTRIBUTORS
%     Alberto Boffi, Enrico Raviola, Andrea Campagna, Luca Ciavirella
% 
% VERSION
%     29-12-2021: v01.0
%-------------------------------------------------------------------------%
    %% INITIAL SETTINGS
    % Position and velocity vectors
    rr = s(1:3);
    vv = s(4:6);

    x= rr(1); y=rr(2); z=rr(3);
    r = norm(rr);
    
    %% ACCELERATIONS DUE TO PERTURBATIONS IN CARTESIAN ELEMENTS
    % gravitational acceleration in cartesian
    a_xyz = [-muP/r^3*rr(1); -muP/r^3*rr(2); -muP/r^3*rr(3)]; % [km/s^2]

    % J2 acceleration in cartesian 
    aJ2_xyz = [  (3/2*J2*muP*Re^2/r^4)*(x/r*(5*z^2/r^2-1));
           (3/2*J2*muP*Re^2/r^4)*(y/r*(5*z^2/r^2-1));
           (3/2*J2*muP*Re^2/r^4)*(z/r*(5*z^2/r^2-3))]; % [km/s^2]

    % SRP acceleration in cartesian
    AU = astroConstants(2); % Astronomical Unit [km]
    pSR = 4.5e-6; % [N/m^2] direct electromagnetic radiation pressure from 
    % the Sun exerted on spacecraft at Earth orbit
    % determination of the Earth's position in heliocentric frame
    [kep, muS] = uplanet(t/(24*60*60), 3);
    % convertion in cartesian elements
    [rE, vE] = kep2car( kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), muS);
    r_sc2s = -(rE + rr); % [km] vector from spacecraft to Sun
    aSRP_xyz = (- pSR * AU^2 * cr * A2m / (norm(r_sc2s))^3 * r_sc2s)/1e+3; % [km/s^2]
  
    % total acceleration
    a_tot_xyz = a_xyz + aJ2_xyz + aSRP_xyz; % [km/s^2]
    

    %% DIFFERENTIAL EQUATIONS OF CARTESIAN ELEMENTS
    ds = [ vv(1); vv(2); vv(3); a_tot_xyz];

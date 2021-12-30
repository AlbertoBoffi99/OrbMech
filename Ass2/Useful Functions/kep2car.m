function [r, v] = kep2car(a, e, i, RAAN, omega, f, muP)
% kep2car Converts keplerian elements of the orbit into the correspondent
% state vector, collecting position and velocity
% 
% 
% PROTOTYPE
%     [r, v] = kep2car(a, e, i, RAAN, omega, f, muP)
%     
% INPUT
%     a [1]                 semi-major axis of the planet     [km]
%     e [1]                 eccentricity magnitude of the planet  [-]
%     i [1]                 inclination of the orbit  [rad]
%     omega [1]             pericenter anomaly of the orbit  [rad]
%     RAAN [1]              right ascension of the ascending node  [-]
%     f [1]                 true anomaly of the body [rad]
%     muP [1]               gravitational parameter of the primary [km^3/s^-2]
%     
% OUTPUT
%     r [3x1]               position vector [km]
%     v [3x1]               velocity vector [km/s]
%
% CONTRIBUTORS
%     Alberto Boffi
% 
% VERSION
%     11-10-2021: v01.0
%
% FUTURE UPGRADES
% add facoltative muP
%-------------------------------------------------------------------------%
    
    % semi-latus rectum
    p = a*(1 - e^2);
    % position vector in peri-focal frame
    rr_pf = [p*cos(f)/(1+e*cos(f)), p*sin(f)/(1+e*cos(f)), 0]';
    % velocity vector in f-rotated peri-focal frame
    vv_rf = [sqrt(muP/p) * e*sin(f), sqrt(muP/p) * (1+e*cos(f)), 0];
    % matrix of rotation of an angle equal to the true anomaly
    Rot_f = [cos(f) sin(f) 0; ...
             -sin(f) cos(f) 0; ...
                 0 0 1];
    % velocity vector in the perifocal frame
    vv_pf = Rot_f' * vv_rf';
    % 3 trivial rotation matrixes
    Rot_RAAN = [cos(RAAN) sin(RAAN) 0; ...
              -sin(RAAN) cos(RAAN) 0; ...
              0 0 1];
    
    Rot_i = [1 0 0; ...
             0 cos(i) sin(i); ...
             0 -sin(i) cos(i)];
    
    Rot_omega = [cos(omega) sin(omega) 0; ...
              -sin(omega) cos(omega) 0; ...
              0 0 1];
    % rotation matrix from Geocentric to peri-focal 
    Rot_GE2PF = Rot_omega * Rot_i * Rot_RAAN;
    % rotation matrix from peri-focal to Geocentric 
    Rot_PF2GE = Rot_GE2PF';
    % position vector in the Geocentric frame
    r = Rot_PF2GE * rr_pf;
    % velocity vector in the Geocentric frame
    v = Rot_PF2GE * vv_pf;
    
end


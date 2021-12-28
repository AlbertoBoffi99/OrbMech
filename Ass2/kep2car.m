function [r,v] = kep2car(a,e,i,OMEGA,omega,theta,muP)
%
% Converter: Keplerian Parameters ---> Cartesian Coordinates & Velocities
%
% DESCRIPTION:
%This code provides a conversion from the Keplerian Parameters of the orbit
%to cartesian coordinates vector (r) and velocities (v).
%--------------------------------------------------------------------------
% INPUTS:
%   a          [1x1]       Semi-Major Axis           [km]
%   e          [1x1]       Eccentricity              [-]
%   i          [1x1]       Inclination               [rad]
%   OM         [1x1]       RAAN                      [rad]
%   om         [1x1]       Argument of Periapsis     [rad]
%   theta      [1x1]       True Anomaly              [rad]
%   muP        [1x1]       Planetary Constant        [km^3][sec-2]
%--------------------------------------------------------------------------
% OUTPUTS:
%   r     [3x1]       Position Vector           [km]
%   v     [3x1]       Velocity Vector           [km/s]
%--------------------------------------------------------------------------

% Rotation matrices: Earth-Centered Inertial --> Perifocal   (ECI->PF)
ROT_OMEGA = [cos(OMEGA) sin(OMEGA) 0;-sin(OMEGA) cos(OMEGA) 0;0 0 1];
ROT_omega = [cos(omega) sin(omega) 0;-sin(omega) cos(omega) 0;0 0 1];
ROT_i = [1 0 0; 0 cos(i) sin(i);0 -sin(i) cos(i)];
%rotation matrix
ROT = ROT_omega*ROT_i*ROT_OMEGA;

% Semi-Latus Rectum
p = a*(1-e^2);

%Position and Velocity Vectors in Perifocal Coordinates
r_pe = [p*(cos(theta)/(1+e*cos(theta))); p*sin(theta)/(1+e*cos(theta)); 0];
v_pe = [-(sqrt(muP/p))*sin(theta);(sqrt(muP/p))*(e+cos(theta));0];

%Position and Velocity Vectors in Cartesian Coordinates
r = ROT'*r_pe;
v = ROT'*v_pe;

end
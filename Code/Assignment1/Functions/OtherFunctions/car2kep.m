function [a, e, i, RAAN, omega, f] = car2kep(r, v, muP)
% car2kep Converts the state vector of the SC into the correspondent keplerian
% elements of the orbit
% 
% PROTOTYPE
%     [a, e, i, RAAN, omega, f] = car2kep(r, v, muP)
%     
% INPUT
%     r [3x1]               position vector [km]
%     v [3x1]               velocity vector [km/s]
%     muP [1]               gravitational parameter of the primary [km^3/s^-2]
%     
% OUTPUT
%     a [1]                 semi-major axis of the planet     [km]
%     e [1]                 eccentricity magnitude of the planet  [-]
%     i [1]                 inclination of the orbit  [rad]
%     omega [1]             pericenter anomaly of the orbit  [rad]
%     RAAN [1]              right ascension of the ascending node  [-]
%     f [1]                 true anomaly of the body [rad]
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

    %% INPUT
    
    % position vector
    rr = r;
    % velocity vector
    vv = v;
    
    %% BODY
    
    % principal versors
    ii = [1,0,0];
    jj = [0,1,0];
    kk = [0,0,1];
    % norm of position vector
    r = norm(rr);
    % norm of velocity vector
    v = norm(vv);
    % vertical velocity
    vr = dot(rr,vv);
    % angular moment vector
    hh = cross(rr,vv);
    % norm of angular velocity
    h = norm(hh);
    % eccentricity vector
    ee = 1/muP .* (cross(vv,hh) - muP.*rr/r);
    % eccentricity norm
    e = norm(ee);
    % semi-major axis
    a = (r * muP)/(2*muP - v^2 *r);
    % node line direction
    NN = cross(kk,hh);
    N = norm(NN);
    % inclination
    i = acos(dot(hh,kk)/h);
    % right ascension of the ascending node
    RAAN = acos(dot(NN,ii)/N);
    if (dot(NN,jj) < 0)
        RAAN = 2*pi - RAAN;
    end
    % pericenter anomaly
    omega = acos(dot(NN,ee)/(N*e));
    if (dot(ee,kk) < 0)
        omega = 2*pi - omega;
    end
    % true anomaly
    f = acos(dot(ee,rr)/(e*r));
    if (vr < 0)
        f = 2*pi - f;
    end
   
    
end

 
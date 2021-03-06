function [dy] = gauss_propagation(t,y,muP,Rp,J2,cr,A2m)
% Propagation of Keplerian elements considering the effect of J2 and SRP 
% 
% Function to be integrated by the ODE solver in order to propagate 
% Keplerian elements in time
% 
% PROTOTYPE
%     [dy] = gauss_propagation(t,y,muP,Rp,J2,cr,A2m)
%     
% INPUT 
%     muP [scalar]   gravitational constant of the planet [km^3/s^2]
%     Rp [scalar]    mean radius of the planet [km]
%     J2 [scalar]    gravitatonal field constant of the Earth
%     cr [scalar]    factor dependent upon the optical properties of
%                      the spacecraft surface
%     A2m[scalar]    surface to mass rate [m^2/kg]
%
% OUTPUT
%     dy [6x1]           Keplerian elements differential equations
%
% CONTRIBUTORS
%     Alberto Boffi, Enrico Raviola, Andrea Campagna, Luca Ciavirella
% 
% VERSION
%     29-12-2021: v01.0
%-------------------------------------------------------------------------%
    
    %% INITIAL SETTINGS
    % Placing the Keplerian elements
    a = y(1); e = y(2); i = y(3); OM = y(4); om = y(5); th = y(6);

    % conversion from Keplerian elements to cartesian
    [r,v] = kep2car(a,e,i,OM,om,th,muP);


    %% ACCELERATIONS DUE TO PERTURBATIONS IN CARTESIAN ELEMENTS
    % J2 related acceleration
    r_norm = norm(r); % [km]
    v_norm = norm(v); % [km/s]
    a_J2_car = ((3*J2*muP*Rp^2)/(2*r_norm^4))* [ (r(1)/r_norm)*(5*(r(3)^2/r_norm^2) - 1); 
                                             (r(2)/r_norm)*(5*(r(3)^2/r_norm^2) - 1); 
                                             (r(3)/r_norm)*(5*(r(3)^2/r_norm^2) - 3)]; % [km/s^2]

    % SRP related acceleration
    psr = 4.5e-6; % [N/m^2] direct electromagnetic radiation pressure from 
    % the Sun exerted on spacecraft at Earth orbit
    AU = astroConstants(2); % Astronomical Unit [km]
    % determination of the Earth's position in heliocentric frame
    [kepE, muS] = uplanet(t/(24*3600), 3);
    % convertion in cartesian elements
    [rE, vE] = kep2car( kepE(1), kepE(2), kepE(3), kepE(4), kepE(5), kepE(6), muS);
    r_sc_S = - (rE + r); % [km] vector from spacecraft to Sun
    a_srp_car =  (- psr * (AU)^2 * cr * A2m / norm(r_sc_S)^3 * r_sc_S)/1e3; % [km/s^2]
    
    a_car = a_J2_car + a_srp_car; % [km/s^2] acceleration only by J2 and SRP
    
    %% DEFINITION OF TNH REFERENCE FRAME
    % Tangential-normal-outofplane RF
    tt = v/norm(v);
    hh = cross(r,v)/norm(cross(r,v));
    nn = cross(hh,tt);
    R_car2tnh = [tt'; nn'; hh']; % Rotation matrix from cartesian to TNH

    %% DETERMINATION OF KEPLERIAN ELEMENTS DIFFERENTIAL EQUATIONS
    a_tnh = R_car2tnh * a_car; % acceleration in TNH frame
    % selecting components
    a_t = a_tnh(1);
    a_n = a_tnh(2);
    a_h = a_tnh(3);

    % Parameters for Gauss planetary equations
    b = a * sqrt(1-e^2);
    n = sqrt(muP/a^3); % orbit's mean motion [rad/s]
    h = n*a*b; % orbit's angular momentum [km^2/s]
    th_star = th + om; % true anomaly + argument of perigee [rad]

    % Gauss planetary equations
    a_dot  = (2*(a^2)*v_norm/muP) * a_t;
    e_dot  = (1/v_norm) * (2*(e+cos(th))*a_t - (r_norm/a)*sin(th)*a_n);
    i_dot  = (r_norm*cos(th_star)/h) * a_h;
    OM_dot = (r_norm*sin(th_star)/(h*sin(i))) * a_h;
    om_dot = (1/(e*v_norm)) * ( 2*sin(th)*a_t + (2*e + (r_norm/a)*cos(th))*a_n ) - (r_norm*sin(th_star)*cos(i)/(h*sin(i)))*a_h;
    th_dot = h/r_norm^2 - (1/(e*v_norm)) * ( 2*sin(th)*a_t + (2*e + (r_norm/a)*cos(th))*a_n );

    % Differential equations of the state (Keplerian elements)
    dy = [a_dot; e_dot; i_dot; OM_dot; om_dot; th_dot];

end
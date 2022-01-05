function p = J2SRP_acc(r,t,cr,A2m);
 
% Earth's gravitational constant
muP = astroConstants(13);
% Earth's radius
Re = astroConstants(23);
% J2 constant
J2 = astroConstants(9);
% radius of the satellite
r_norm = sqrt(dot(r,r));
% SRP related acceleration
psr = 4.5e-6; % [N/m^2] direct electromagnetic radiation pressure from 
% the Sun exerted on spacecraft at Earth orbit
AU = astroConstants(2); % Astronomical Unit [km]
    
    % acceleration from secular disturbance
    s = [(3/2*J2*muP*Re^2/r_norm^4)*(r(1)/r_norm*(5*r(3)^2/r_norm^2-1));
       (3/2*J2*muP*Re^2/r_norm^4)*(r(2)/r_norm*(5*r(3)^2/r_norm^2-1));
       (3/2*J2*muP*Re^2/r_norm^4)*(r(3)/r_norm*(5*r(3)^2/r_norm^2-3))];

    % determination of the Earth's position in heliocentric frame
    [kepE, muS] = uplanet(t/(24*3600), 3);
    % convertion in cartesian elements
    [rE, ~] = kep2car( kepE(1), kepE(2), kepE(3), kepE(4), kepE(5), kepE(6), muS);
    % distance between SC and sun
    r_S = - (rE + r); 
    q =  (- psr * (AU)^2 * cr * A2m / norm(r_S)^3 * r_S)/1e3; 

    p = s + q;

end

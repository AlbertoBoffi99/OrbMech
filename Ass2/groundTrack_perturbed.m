% function to calculate groundtracks considering J2 perturbation

function [alpha, delta, lon, lat] = groundTrack_perturbed (K_elements, lon0, tspan, omeP, muP, Re, J2, cr, A2m, t0)

%% function 

[r0,v0] = kep2car(K_elements(1), K_elements(2), K_elements(3), K_elements(4), K_elements(5), K_elements(6), muP);

C_elements = [r0; v0];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

[ ~, state ] = ode113( @(t,s) ode_perturbed(t,s,muP,Re,J2,cr,A2m), tspan, C_elements, options );

r = sqrt((state(:,1)).^2 + (state(:,2)).^2 + (state(:,3)).^2);

% declination in rad
delta = asin(state(:,3)./r);

% right ascension in rad
alpha = [];
for k = 1:length(tspan)
    
    if state(k,2)/r(k) > 0
        alpha(k) = acos((state(k,1)/r(k))/cos(delta(k)));
    else
        alpha(k) = 2*pi - acos((state(k,1)/r(k))/cos(delta(k)));
    end

end 

lon_E = lon0 + omeP*(tspan-t0);

lon = alpha - lon_E;

lat = delta;

end
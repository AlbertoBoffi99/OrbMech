% function to calculate a (semimajor axis) from k and m given integers
% PERTURBED

function a_r_p = repeatingGroundTracks_perturbed(k, m, muP, omeP, Re, J2, e, i)

% % longitude of the ascending Node rate
% Ome_dot = @(a) -((3*sqrt(muP)*J2*Re^2)/(2*(1-e^2)^2*a^(7/2)))*cos(i);
% 
% % argument of pheriapsis rate 
% ome_dot = @(a) -((3*sqrt(muP)*J2*Re^2)/(2*(1-e^2)^2*a^(7/2)))*(5/2*(sin(i))^2-2);
% 
% % mean anomaly rate
% M_dot = @(a) ((3*sqrt(muP)*J2*Re^2)/(2*((1-e^2)^(3/2))*a^(7/2)))*(1-3/2*(sin(i))^2);
% 
% % angular rate in a circular orbit (constant velocity)
% n = @(a) sqrt(muP/a^3);

% a_solve to obtain roots
a_solve = @(a) (m/k - (omeP - (-((3*sqrt(muP)*J2*Re^2)/(2*(1-e^2)^2*a^(7/2)))*cos(i)))/(sqrt(muP/a^3) - ((3*sqrt(muP)*J2*Re^2)/(2*(1-e^2)^2*a^(7/2)))*(5/2*(sin(i))^2-2) + ((3*sqrt(muP)*J2*Re^2)/(2*(1-e^2)^(3/2)*a^(7/2)))*(1-3/2*(sin(i))^2)));

f_zero_options = optimset('TolFun',1e-12);

% a_guess
a_guess = repeatingGroundTracks(k,m,muP,omeP);

% a caught from fsolve function
a_r_p = fzero( a_solve , a_guess, f_zero_options );


end
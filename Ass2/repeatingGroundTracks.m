function a_r = repeatingGroundTracks(k, m, muP, omeP)
% Function to calculate a (semimajor axis) of the satellite's unperturbed
% orbit that leads to a repeating ground track
% 
% PROTOTYPE
%     a_r = repeatingGroundTracks(k, m, muP, omeP)
%     
% INPUT 
%     k [scalar]        number of satellite's revolutions 
%     m [scalar]        number of planet's rotations
%     muP [scalar]      gravitational constant of the planet [km^3/s^2]
%     omeP [scalar]     planet's rotational angular velocity [rad/s]
%
% OUTPUT
%     a_r [scalar]      repeating ground track orbit's semimajor axis [km]
%
% CONTRIBUTORS
%     Alberto Boffi, Enrico Raviola, Andrea Campagna, Luca Ciavirella
% 
% VERSION
%     29-12-2021: v01.0
%-------------------------------------------------------------------------%
%% INITIAL SETTINGS
% satellite angular velocity
n = k/m*omeP; % [rad/s]

%% OUTPUT
a_r = (muP/n^2)^(1/3);

end
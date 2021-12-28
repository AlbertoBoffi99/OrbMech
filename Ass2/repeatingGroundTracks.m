% function to calculate a (semimajor axis) from k and m given integers

function a_r = repeatingGroundTracks(k, m, muP, omeP)

% satellite angular velocity
n = k/m*omeP; % [rad/s]

% semimajor axis for specific inputs
a_r = (muP/n^2)^(1/3); % [rad]

end
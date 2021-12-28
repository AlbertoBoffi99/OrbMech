% Differential equations for the two-body problem

% Function to compute the motion of a secondary object around a primary one
% with much bigger mass, under the only effect of the gravitaional force
% and no other influences

function dy = tbp_ode ( ~, y, muP ) 

%% INPUT

% position components
r = y(1:3); % [km]

% velocity components
v = y(4:6); % [km/s]

% distance of m1 from m2
r_norm = norm(r); % [km]

%% BODY

% system of first oreder ODEs which gives the derivative of 
% the state

dy = [ v(1); v(2); v(3); 
    -muP/r_norm^3*r(1)
    -muP/r_norm^3*r(2)
    -muP/r_norm^3*r(3)
    ];

end
    


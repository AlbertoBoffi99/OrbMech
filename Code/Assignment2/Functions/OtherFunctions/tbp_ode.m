function ds = tbp_ode(~, s, muP)
% tbp_ode Differential equations for the two-body problem
% 
% Function to compute the motion of a secondary object around a primary one
% with much bigger mass, under the only effect of the grivitaional force
% and no other influences
% 
% PROTOTYPE
%     ds = tbp_ode(~, s, muP)
%     
% INPUT
%     ~ [1]           time [s] (SUPPRESSED)
%     s [6x1]         state vector [rx ry rz vx vy vz]'
%                       where
%                       ri: i-component of position     [km]
%                       vi: i-component of velocity     [km/s]
%     muP [1]         gravitational parameter of the primary [km^3/s^-2]
%     
% OUTPUT
%     ds [1x6]        differential of state vector [vx vy vz ax ay az]'   
%                       where
%                       vi: i-component of velocity     [km/s]
%                       ai: i-component of acceleration [km/s^2]
% CONTRIBUTORS
%     Alberto Boffi
% 
% VERSION
%     24-09-2021: v01.0
%-------------------------------------------------------------------------%

%% INPUT

    % position components
    r = s(1:3); %[km]
    % velocity components
    v = s(4:6); %[km/s]
    % distance from m1
    r_norm = norm(r);   %[km]
    
%% BODY
    
    % system of first oreder ODEs which gives the derivative of the state
    % vector s
    ds = [
        v(1)
        v(2)
        v(3)
        -muP/r_norm^3*r(1)
        -muP/r_norm^3*r(2)
        -muP/r_norm^3*r(3)
        ];
    
end
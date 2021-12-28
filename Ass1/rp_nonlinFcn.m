function [c, ceq] = rp_nonlinFcn(x, RE, h_atm)
% rp_nonlinFcn Non linear constraints for GA or fmincon optimizers
% 
% Function to create a non linear constraint on the GA or fmincon
% optimization which checks whether the perigee radius is greater than the
% sum of atmpsphere altitude and Earth's radius or not
% 
% PROTOTYPE
%     [c, ceq] = rp_nonlinFcn(x, RE, h_atm)
%     
% INPUT
%     x [nx1]           vector of optimization variables used also in
%                       optimizer call
%     RE [1]            Earth's radius
%     h_atm [1]         atmosphere height
%     
% OUTPUT
%     c [mx1]           array of nonlinear inequality constraints at x
%     ceq [fx1]         array of nonlinear equality constraints at x
%
% CONTRIBUTORS
%     Alberto Boffi, Enrico Raviola, Andrea Campagna, Luca Ciavirella
% 
% VERSION
%     26-12-2021: v01.0
%-------------------------------------------------------------------------%

    %% INPUT

    % calling global variable
    global out

    %% INEQUALITY CONSTRAINTS

    % perigee radius has to be greater than Earth's radius + atmosphere's
    % height

    if out.rp_norm < RE + h_atm
        % if this is not true the constraint value is as high as the
        % difference between the perigee radius and the minimum radius to
        % satisfy
        c = abs(out.rp_norm - (RE + h_atm));
    else
        c = -1;
    end

     %% EQUALITY CONSTRAINTS
    
     % no non-linear inequality constraint is set
     ceq = [];

end
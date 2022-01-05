function [A] = DCM312 (theta, phi, psi)
% DCM312 Evaluates the DCM with the 312 covention using Euler angles 
% 
% PROTOTYPE
%     [A] = DCM312 (theta, phi, psi)
%     
% INPUT
%     theta [1]             theta Euler angle
%     phi [1]               phi Euler angle
%     psi [1]               psi Euler angle
%
%     
% OUTPUT
%     A [3x3]               Direction Cosine Matrix
%
% CONTRIBUTORS
%     Alberto Boffi
% 
% VERSION
%     20-11-2021: v01.0
%
%-------------------------------------------------------------------------%

A = [cos(psi)*cos(phi)-sin(psi)*sin(phi)*sin(theta) cos(psi)*sin(phi)+sin(psi)*cos(phi)*sin(theta) -sin(psi)*cos(theta);...
           -sin(phi)*cos(theta)                           cos(phi)*cos(theta)                            sin(theta);...
           sin(psi)*cos(phi)+cos(psi)*sin(phi)*sin(theta) sin(psi)*sin(phi)-cos(psi)*cos(phi)*sin(theta) cos(theta)*cos(psi)];
       
end
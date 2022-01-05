function [A] = DCM313 (theta, phi, psi)
% DCM312 Evaluates the DCM with the 312 covention using Euler angles 
% 
% PROTOTYPE
%     [A] = DCM313 (theta, phi, psi)
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

A = [cos(psi)*cos(phi)-sin(psi)*sin(phi)*cos(theta)  cos(psi)*sin(phi)+sin(psi)*cos(phi)*cos(theta)  sin(psi)*sin(theta);...
           -sin(psi)*cos(phi)-cos(psi)*sin(phi)*cos(theta) -sin(psi)*sin(phi)+cos(psi)*cos(phi)*cos(theta) sin(theta)*cos(psi);...       
           sin(phi)*sin(theta)                             -cos(phi)*sin(theta)                            cos(theta)];
       
end
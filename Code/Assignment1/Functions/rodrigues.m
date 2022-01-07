function [V] = rodrigues(v, u, delta)
% rodrigues Rodrigues vector rotation
% 
% Function to compute rotation of a vector through Rodrigues method
% 
% PROTOTYPE
%     [v_plus] = rodrigues(v, u, delta)
%     
% INPUT
%     v [3x1]       vector to rotate
%     u [3x1]       vector around which v vector is rotated
%     delta [1]     angle of which v vector is rotated
%     
% OUTPUT
%     V [3x1]         rotated v vector
%
% CONTRIBUTORS
%     Alberto Boffi, ...
% 
% VERSION
%     16-12-2021: v01.0
%-------------------------------------------------------------------------%

    V = v * cos(delta) + cross(u,v) * sin(delta)+ u * dot(u,v) * (1-cos(delta));

end
function [rp_norm, rp, DvGA, Dvfb, vp_p, vp_m] = GAflyby(V_min, V_plus, v_pl, mu_pl, R_pl, h_atm_pl, options)
% GAflyby Assisted flyby computation
% 
% Function to compute the assisted flyby around a planet
% 
% PROTOTYPE
%     [rp_norm, rp, Dv_p, vp_p, vp_m] = GAflyby(V_min, V_plus, v_pl, mu_pl, R_pl, h_atm_pl, options)
%     
% INPUT
%     V_min [3x1]       absolute velocity vector at minus infinte
%     V_plus [3x1]      absolute velocity vector at plus infinte
%     v_pl [3x1]        planet velocity vector at flyby instant
%     mu_pl [1]         planet gravitational constant
%     R_pl [1]          fly-by planet radius
%     h_atm_pl [1]      planet's atmosphere height
%     options [struct]  fsolve options
%     
% OUTPUT
%     rp_norm [1]       perigee radius vecto norm
%     rp [3x1]          perigee radius vector
%     DvGA[1]           delta velocity given at perigee
%     Dvfb[1]           natural fly-by delta velocity
%     vp_p[3x1]         perigee velocity vector on exit hyperbola         
%     vp_m[3x1]         perigee velocity vector on entry hyperbola 
%
% CONTRIBUTORS
%     Alberto Boffi, Enrico Raviola, Andrea Campagna, Luca Ciavirella
% 
% VERSION
%     16-12-2021: v01.0
%-------------------------------------------------------------------------%


    %% FLY-BY CALCULATION
    
    % infinity relative velocity at entry
    v_min = [V_min(1)-v_pl(1) V_min(2)-v_pl(2) V_min(3)-v_pl(3)];
    % infinity relative velocity ar exit
    v_plus = [V_plus(1)-v_pl(1) V_plus(2)-v_pl(2) V_plus(3)-v_pl(3)];

    % natural fly-by delta velocity
    Dvfb = norm(v_plus - v_min);

    % turning angle
    TA = acos(dot(v_min,v_plus)/(norm(v_min)*norm(v_plus)));

    % eccentricity and turning angle of entry hyperbola
    e_min = @(rp_norm) 1+ (rp_norm*(norm(v_min))^2)/mu_pl;
    TA_min = @(rp_norm) 2*asin(1/e_min(rp_norm));

    % eccentricity and turning angle of exit hyperbola
    e_plus = @(rp_norm) 1+ (rp_norm*(norm(v_plus))^2)/mu_pl;
    TA_plus = @(rp_norm) 2*asin(1/e_plus(rp_norm));

    % implicit function to infer perigee radius
    rp_fun = @(rp) TA_min(rp)+TA_plus(rp)-2*TA;

    % perigee radius
    rp_norm = fsolve(rp_fun, R_pl + h_atm_pl, options);  

    % eccentricity and turning angle of entry hyperbola
    e_min = 1+ (rp_norm*(norm(v_min))^2)/mu_pl;
    TA_min = 2*asin(1/e_min);

    % eccentricity and turning angle of exit hyperbola
    e_plus = 1+ (rp_norm*(norm(v_plus))^2)/mu_pl;
    TA_plus = 2*asin(1/e_plus);

    % perigee relative velocity at entry velocity
    vp_min = sqrt((norm(v_min)^2) + 2*mu_pl/rp_norm);
    % perigee relative velocity at exit velocity
    vp_plus = sqrt((norm(v_plus)^2) + 2*mu_pl/rp_norm);                
    
    % rotation vector for rodrigues function
    u = cross(v_min,v_plus)/norm(cross(v_min,v_plus));
    
    % perigee velocity direction on entry hyperbola
    dir_vp_min = rodrigues(v_min, u, TA_min/2)/norm(v_min);
    % perigee velocity direction on exit hyperbola
    dir_vp_plus = rodrigues(v_plus, -u, TA_plus/2)/norm(v_plus);  
    
    % perigee velocity vector on entry hyperbola
    vp_m = vp_min * dir_vp_min;      
    % perigee velocity vector on exit hyperbola
    vp_p = vp_plus * dir_vp_plus;      
    
    % perigee radius direction
    dir_rp = cross( vp_m, u)/norm(cross(vp_m,u));  
    % perigee radius vector
    rp = rp_norm * dir_rp;
    
    % flyby delta velocity at perigee
    DvGA = norm(vp_p - vp_m);


end




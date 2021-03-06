function [Dv_tot] = JEV_lamb_fb_lamb(dates, astro, options)
% JEV_lamb_fb_lamb Jupiter - Venus with one Earth assisted fly-by
% 
% Function to compute the interplanetary trajectory:
% 1. Lambert arc Jupiter - Earth
% 2. Earth assisted fly-by
% 3. Lambert arc Earth - Venus
% 
% PROTOTYPE
%     [Dv_tot] = JEV_lamb_fb_lamb(dates, astro, options)
%     
% INPUT
%     dates [3]         vector of mjd2000 dates: departure, flyby, arrival
%     astro [struct]    structure containing astronomical constant used in
%                       INPUT section
%     options [struct]  structure containing Matlab function options used in
%                       INPUT section
%     
% OUTPUT
%     Dv_tot [3x1]      total delta velocity vector
%
% CONTRIBUTORS
%     Alberto Boffi, Enrico Raviola, Andrea Campagna, Luca Ciavirella
% 
% VERSION
%     16-12-2021: v01.0
%-------------------------------------------------------------------------%
    
    %% INPUT

    % calling global function out
    global out 

    % astro
    muS = astro.muS;
    muE = astro.muE;
    RE = astro.RE;
    h_atm = astro.h_atm;
    % options
    orbitType = options.orbitType;
    Nrev = options.Nrev;
    Ncase = options.Ncase;
    LambOptions = options.LambOptions;
    fsolveOptions = options.fsolve_options;
    
    %% INITIAL PLANET STATE
    
    % Jupiter departure state
    kepJ = uplanet(dates(1), 5);
    [rdJ, vdJ]= kep2car(kepJ(1),kepJ(2),kepJ(3),kepJ(4),kepJ(5),kepJ(6),muS);
    % Earth flyby state
    kepE = uplanet(dates(2), 3);
    [rfbE, vfbE]= kep2car(kepE(1),kepE(2),kepE(3),kepE(4),kepE(5),kepE(6),muS);
    % Venus arrival state
    kepV= uplanet(dates(3), 2);
    [raV, vaV]= kep2car(kepV(1),kepV(2),kepV(3),kepV(4),kepV(5),kepV(6),muS);

    %% FIRST LAMBERT ARC
    
    % time of flight of 1st arc
    T1 = (dates(2) - dates(1))*24*60*60;
 
    % lambert problem for 1st arc
    [~,~,~,~,ViT1,VfT1,~,~] = lambertMR(rdJ,rfbE,T1,muS,orbitType,Nrev,Ncase,LambOptions);
    
    % 1st arc delta velocity
    DviT1x= ViT1(1)-vdJ(1);
    DviT1y= ViT1(2)-vdJ(2);
    DviT1z= ViT1(3)-vdJ(3);
    % norm of 1st arc delta velocity
    DviT1= norm([DviT1x DviT1y DviT1z]);
    % renaming delta velocity of 1st arc
    DvT1 = DviT1;

    %% SECOND LAMBERT ARC

    % time of flight of 2nd arc
    T2 = (dates(3) - dates(2))*24*60*60;

    % lambert problem for 2nd arc
    [~,~,~,~,ViT2,VfT2,~,~] = lambertMR(rfbE,raV,T2,muS,orbitType,Nrev,Ncase,LambOptions);
    
    % 2nd arc delta velocity
    DvfT2x = vaV(1)-VfT2(1);
    DvfT2y = vaV(2)-VfT2(2);
    DvfT2z = vaV(3)-VfT2(3);
    % norm of 2nd arc delta velocity
    DvfT2 = norm([DvfT2x DvfT2y DvfT2z]);
    % renaming delta velocity of 2nd arc
    DvT2 = DvfT2;

    %% FLY-BY
    
    % performing fly-by
    [rp_norm, rp, DvGA, Dvfb, vp_p, vp_m] = GAflyby(VfT1, ViT2, vfbE, muE, RE, h_atm, fsolveOptions);

    % total delta velocity
    Dv_tot = DvT1 + DvT2 + DvGA;
    out.nanflag = 0;

    % check on perigee radius
    % if the perigee radius is higher than the atmosphere of the planet
    % then save all useful  outputs, otherwise save only discarded values
    % and a flag variable
    if rp_norm > RE + h_atm
        out.ViT1 = ViT1;
        out.ViT2 = ViT2;
        out.DvT1 = DvT1;
        out.DvT2 = DvT2;
        out.DvGA = DvGA;
        out.Dvfb = Dvfb;
        out.rp_norm = rp_norm;
        out.rp = rp;
        out.vp_p = vp_p;
        out.vp_m = vp_m;
    else
        out.Dv_disc = DvT1 + DvT2 + DvGA;
        out.dates_disc = dates;
        out.nanflag = 1;
        out.rp_norm = rp_norm;
        Dv_tot = 1e4;
        % warning("Perigee radius is below Earth's atmosphere")
    end


end
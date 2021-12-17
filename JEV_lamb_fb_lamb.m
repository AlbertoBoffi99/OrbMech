function Dv_tot = JEV_lamb_fb_lamb(dates, astro, options)
% JEV_lamb_fb_lamb Jupiter - Venus with one Earth assisted fly-by
% 
% Function to compute the interplanetary trajectory:
% 1. Lambert arc Jupiter - Earth
% 2. Earth assisted fly-by
% 3. Lambert arc Earth - Venus
% 
% PROTOTYPE
%     Dv_tot = JEV_lamb_fb_lamb(dates, astro, options)
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
%     Alberto Boffi, ...
% 
% VERSION
%     16-12-2021: v01.0
%-------------------------------------------------------------------------%
    
    %% INPUT

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
    
    % Jupiter initial state
    kepJ = uplanet(dates(1), 5);
    [rdJ, vdJ]= kep2car(kepJ(1),kepJ(2),kepJ(3),kepJ(4),kepJ(5),kepJ(6),muS);
    % Earth initial state
    kepE = uplanet(dates(2), 3);
    [rfbE, vfbE]= kep2car(kepE(1),kepE(2),kepE(3),kepE(4),kepE(5),kepE(6),muS);
    % Venus initial state
    kepV= uplanet(dates(3), 2);
    [raV, vaV]= kep2car(kepV(1),kepV(2),kepV(3),kepV(4),kepV(5),kepV(6),muS);

    %% FIRST LAMBERT ARC
    
    % time of flight of 1st arc
    T1 = (dates(2) - dates(1))*24*60*60;  
    % lambert problem for 1st arc
    [aT1,pT,E_T,ERROR,ViT1,VfT1,TPAR,THETA] = lambertMR(rdJ,rfbE,T1,muS,orbitType,Nrev,Ncase,LambOptions);
    
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
    [aT2,pT,E_T,ERROR,ViT2,VfT2,TPAR,THETA] = lambertMR(rfbE,raV,T2,muS,orbitType,Nrev,Ncase,LambOptions);
    
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
    [rp_norm, rp, Dvfb] = GAflyby(VfT1, ViT2, vfbE, muE, RE, h_atm, fsolveOptions);
    
    % total delta velocity
    Dv_tot= DvT1 + DvT2 + Dvfb ;


end
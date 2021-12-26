%% PLOT INTERPLANETARY TRAJECTORY

% Jupiter's position at selected departure data
kepJ = uplanet(results.dates(1), 5);
[rdJ, vdJ]= kep2car(kepJ(1),kepJ(2),kepJ(3),kepJ(4),kepJ(5),kepJ(6),astro.muS);
% Earth's position at selected flyby date
kepE = uplanet(results.dates(2), 3);
[rfbE, vfbE]= kep2car(kepE(1),kepE(2),kepE(3),kepE(4),kepE(5),kepE(6),astro.muS);
% Venus' position at selected arrival date
kepV = uplanet(results.dates(3), 2);
[raV, vaV]= kep2car(kepV(1),kepV(2),kepV(3),kepV(4),kepV(5),kepV(6),astro.muS);
% Jupiter's orbit propagation
TJ = 2*pi*sqrt( (kepJ(1))^3/astro.muS );              
tspanJ = linspace( 0, TJ, 1000 );
[~, stateJ] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspanJ, [rdJ, vdJ], options.ode_options);
% Earth's orbit propagation
TE = 2*pi*sqrt( (kepE(1))^3/astro.muS );             
tspanE = linspace( 0, TE, 1000 );
[~, stateE] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspanE, [rfbE, vfbE], options.ode_options);
% Venus' orbit propagation
TV = 2*pi*sqrt( (kepV(1))^3/astro.muS );              
tspanV = linspace( 0, TV, 1000 );
[~, stateV] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspanV, [raV, vaV], options.ode_options);

% Lambert arc 1 and 2 propagation
T1 = (results.dates(2)-results.dates(1))*24*60*60;
T2 = (results.dates(3)-results.dates(2))*24*60*60;

tspan_t1 = linspace( 0, T1, 1000);
[~, state_t1] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspan_t1, [rdJ, results.ViT1'], options.ode_options);
tspan_t2 = linspace( 0, T2, 1000 );
[~, state_t2] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspan_t2, [rfbE, results.ViT2'], options.ode_options);

figure()  
plot3( stateJ(:,1), stateJ(:,2), stateJ(:,3), 'red', 'linewidth', 2)        
hold on
plot3( state_t1(:,1), state_t1(:,2), state_t1(:,3), 'green--', 'linewidth', 1)
plot3( stateE(:,1), stateE(:,2), stateE(:,3), 'blue', 'linewidth', 2)     
plot3( state_t2(:,1), state_t2(:,2), state_t2(:,3), 'magenta--', 'linewidth', 1) 
plot3( stateV(:,1), stateV(:,2), stateV(:,3), 'yellow', 'linewidth', 2)     
% add Planets at corresponding positions
opts.units= 's';
planet3D('Sun', opts, [0 0 0]);  
opts.units= 'j';
planet3D( 'Jupiter', opts, rdJ);
opts.units= 'p';
planet3D('Earth', opts, rfbE);
planet3D( 'Venus', opts, raV);

legend('Jupiter orbit','Transfer arc 1', 'Earth orbit','Transfer arc 2','Venus orbit', 'Sun', 'Jupiter at launch','Earth at flyby', 'Venus at arrival');

xlabel("r_x[km]");
ylabel("r_y[km]");
zlabel("r_z[km]");
title('Trajectory');
grid on
axis equal

%% PRELIMINARY ANALYSIS OF FIRST LAMBERT ARC

in_dep_date= date2mjd2000([2025 08 01 00 00 00]);  
fin_dep_date= date2mjd2000([2060 08 01 00 00 00]);  

in_fb_date= date2mjd2000([ 2025 08 01 00 00 00]);
fin_fb_date= date2mjd2000([ 2061 08 01 00 00 00]);

n=1000;

dep_window= linspace(in_dep_date, fin_dep_date,n);
fb_window= linspace(in_fb_date, fin_fb_date,n);

Dvmin = astroConstants(5);                                
r = NaN;
c = NaN;                                                  
% select an interval to cut the contour plot
piano = (1:1:5);                                    

for i = 1 : n
    % for every Jupiter's position at departure
    kep1 = uplanet(dep_window(i), 5);
    [r1(:,i), v1(:,i)] = kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),astro.muS);

    for j = 1 : n
        %for every Earth's position at flyby
        kep2 = uplanet(fb_window(j), 3);
        [r2(:,j), v2(:,j)] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),astro.muS);
        
        T11(i,j) = (fb_window(j)-dep_window(i))*24*60*60;   
        %flyby corresponding to the i-departure date and j-flyby date
        [aT1(i,j),pT,E_T,ERROR,Vi1(i,j,:),Vf1(i,j,:),TPAR,THETA] = ...
            lambertMR(r1(:,i),r2(:,j),T11(i,j),astro.muS,options.orbitType,options.Nrev,options.Ncase,options.LambOptions);
        Dv1x = Vi1(i,j,1)-v1(1,i);
        Dv1y = Vi1(i,j,2)-v1(2,i);
        Dv1z = Vi1(i,j,3)-v1(3,i);
        Dv1 = norm([Dv1x Dv1y Dv1z]); %delta velocity to insert into arc 1   
        Dvtot(i,j)= Dv1;    
        %we don't need to minimize delta velovity at Earth's arrival since
        %we're perfoming a flyby
        
        if Dvtot(i,j)<Dvmin 
            Dvmin= Dvtot(i,j);
            r=i;
            c=j;
        end
    end
    x(i)=datenum(mjd20002date(dep_window(i)));  
    y(i)=datenum(mjd20002date(fb_window(i))); 
end 
% optimum dates for single arc
launch_day= mjd20002date(dep_window(r));
flyby_day= mjd20002date(fb_window(c));
% dates obtained from overall problem
l_day= mjd20002date(results.dates(1));
fb_day= mjd20002date(results.dates(2));

figure()
hold on
[C, h] = contour(x, y, Dvtot', floor(Dvmin) + (piano));
plot(x(r),y(c), '*', 'LineWidth',20);
caxis([Dvmin Dvmin + piano(length(piano))]);                                     
Color= colorbar;                                              
clabel(C,h, floor(Dvmin) + (piano));                         
datetick( 'x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
datetick( 'y', 'yyyy mmm dd', 'keeplimits');
xtickangle(45);
xlabel('Jupiter departure date');
ylabel('Earth flyby date');
set(Color,'YtickLabel');
Color.Label.String= 'Delta V1 [km/s]';
grid on
title('Pork Chop Plot');

plot(datenum(l_day),datenum(fb_day), 'd', 'LineWidth',20);
scatter(datenum(mjd20002date(results.dates(1))), datenum(mjd20002date(results.dates(2))), 'LineWidth', 12);
% plot time of flight in the same contour
[C1, h1] = contour(x, y, (T11/(24*3600))','Black');
clabel(C1,h1);                                              

%% PRELIMINARY ANALYSIS OF SECOND LAMBERT ARC

in_arr_date= date2mjd2000([2027 08 01 00 00 00]);
fin_arr_date= date2mjd2000([2065 08 01 00 00 00]);  

n=1000;      
arr_window= linspace(in_arr_date, fin_arr_date,n);

Dvmin=astroConstants(5);                             
r=NaN;
c=NaN;                                                  

piano = (1:1:10);                                

for j=1:n
    kep2= uplanet(fb_window(j), 3);
    [r2(:,j), v2(:,j)]= kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),astro.muS);

    for k=1:n
        kep3= uplanet(arr_window(k), 2);
        [r3(:,k), v3(:,k)]= kep2car(kep3(1),kep3(2),kep3(3),kep3(4),kep3(5),kep3(6),astro.muS);
        
        T22(j,k)= (arr_window(k)-fb_window(j))*24*60*60;   

        [aT2(j,k),pT,E_T,ERROR,Vi2(j,k,:),Vf2(j,k,:),TPAR,THETA] = ...
            lambertMR(r2(:,j),r3(:,k),T22(j,k),astro.muS,options.orbitType,options.Nrev,options.Ncase,options.LambOptions);
        Dv2x= v3(1,k)-Vf2(j,k,1);
        Dv2y= v3(2,k)-Vf2(j,k,2);
        Dv2z= v3(3,k)-Vf2(j,k,3);
        Dv2= norm([Dv2x Dv2y Dv2z]);   
        
        Dvtot(j,k)= Dv2;    

        if Dvtot(j,k)<Dvmin 
            Dvmin= Dvtot(j,k);
            r=j;
            c=k;
        end
    end
    x(j) = datenum(mjd20002date(fb_window(j)));  
    y(j) = datenum(mjd20002date(arr_window(j))); 
end 
% optimum dates for single arc
flyby_day= mjd20002date(fb_window(r));
arrival_day= mjd20002date(arr_window(c));
% dates obtained from overall problem
fb_day= mjd20002date(results.dates(2));
a_day= mjd20002date(results.dates(3));

figure()
hold on
[C, h] = contour(x, y, Dvtot', floor(Dvmin) + (piano));
plot(x(r),y(c), '*', 'LineWidth',20);
caxis([Dvmin Dvmin + piano(length(piano))]);                                     
Color= colorbar;                                              
clabel(C,h, floor(Dvmin) + (piano));                         
datetick( 'x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
datetick( 'y', 'yyyy mmm dd', 'keeplimits');
xtickangle(45);
xlabel('Earth flyby date');
ylabel('Venus arrival date');
set(Color,'YtickLabel');
Color.Label.String= 'Delta V2 [km/s]';
grid on
title('Pork Chop Plot');
plot(datenum(fb_day),datenum(a_day), 'd', 'LineWidth',20);
scatter(datenum(mjd20002date(results.dates(2))), datenum(mjd20002date(results.dates(3))), 'LineWidth', 12);


[C1, h1] = contour(x, y, (T22/(24*3600))','Black');
clabel(C1,h1);                                              

%% EARTH'S FLYBY

tspanH = linspace( 0, 5000, 1000 );
% outgoing hyperpola propagation
[timeH1, stateH1] = ode113( @(t,s) tbp_ode(t, s, astro.muE), tspanH, [results.rp,results.vp_p], options.ode_options);
%incoming hyperbola propagation
[timeH2, stateH2] = ode113( @(t,s) tbp_ode(t, s, astro.muE), -tspanH, [results.rp,results.vp_m], options.ode_options); 
r_H1= stateH1(:,1:3);
r_H2= stateH2(:,1:3);
%possono essere aggiunte in config
G = astroConstants(1);
ME = astroConstants(13)/G;
MS = astroConstants(4)/G;

r_soi = astroConstants(2)*(ME/MS)^(2/5);
flag1=1;
flag2=1;

for kk=1:length(tspanH)
    if norm(r_H1(kk,:)) > r_soi && flag1==1
        time1 = timeH1(kk);
        flag1 = 0;
    end
    if norm(r_H2(kk,:)) > r_soi && flag2==1
        time2 = -timeH2(kk);
        flag2 = 0;
    end
end
%time spent inside the Earth's SOI
fb_time = time1 + time2;

figure()  
p1 = plot3( stateH1(:,1), stateH1(:,2), stateH1(:,3), 'red', 'linewidth', 2);   %outgoing hyperbola  
hold on
p2 = plot3( stateH2(:,1), stateH2(:,2), stateH2(:,3), 'blue', 'linewidth', 2);  %incoming hyperbola
opts.units = 'km';
planet3D('Earth',opts, [0 0 0]);    %add planet
grid on
axis equal
% add arrow in the direction of the planet's velocity 
p3 = arrow3d([0 vfbE(1)*10^3],[0 vfbE(2)*10^3],[0 vfbE(3)*10^3],9/10,0.2*10^3,0.8*10^3,'black'); 
% showing the perigee passage
p4 = scatter3(stateH1(1,1), stateH1(1,2),stateH1(1,3), 'magenta' ,'filled', 'LineWidth', 10);

title('Flyby in Earth-centred frame parallel to HECI')
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
legend([p2 p1 p3 p4],'incoming hyperbola', 'outgoing hyperbola', 'planet velocity', 'pericenter' )

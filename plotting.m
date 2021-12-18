muS = astro.muS; 
muE = astro.muE;

kepJ = uplanet(results.fmincon_dates(1), 5);
[rdJ, vdJ]= kep2car(kepJ(1),kepJ(2),kepJ(3),kepJ(4),kepJ(5),kepJ(6),muS);

kepE = uplanet(results.fmincon_dates(2), 3);
[rfbE, vfbE]= kep2car(kepE(1),kepE(2),kepE(3),kepE(4),kepE(5),kepE(6),muS);

kepV = uplanet(results.fmincon_dates(3), 2);
[raV, vaV]= kep2car(kepV(1),kepV(2),kepV(3),kepV(4),kepV(5),kepV(6),muS);

options = odeset ('RelTol', 1e-3, 'AbsTol', 1e-3);

TJ = 2*pi*sqrt( (kepJ(1))^3/muS );              %period of Earth around the Sun
tspanJ = linspace( 0, TJ, 1000 );
[timeJ, stateJ] = ode113( @(t,s) tbp_ode(t, s, muS), tspanJ, [rdJ, vdJ], options);

TE = 2*pi*sqrt( (kepE(1))^3/muS );              %period of Earth around the Sun
tspanE = linspace( 0, TE, 1000 );
[timeE, stateE] = ode113( @(t,s) tbp_ode(t, s, muS), tspanE, [rfbE, vfbE], options);

TV = 2*pi*sqrt( (kepV(1))^3/muS );              %period of Earth around the Sun
tspanV = linspace( 0, TV, 1000 );
[timeV, stateV] = ode113( @(t,s) tbp_ode(t, s, muS), tspanV, [raV, vaV], options);


T1= (results.fmincon_dates(2)-results.fmincon_dates(1))*24*60*60;
T2= (results.fmincon_dates(3)-results.fmincon_dates(1))*24*60*60;

[aT1,pT,E_T,ERROR,ViT1,VfT1,TPAR,THETA] = lambertMR(rdJ,rfbE,T1,muS,orbitType,Nrev,Ncase,LambOptions);
[aT2,pT,E_T,ERROR,ViT2,VfT2,TPAR,THETA] = lambertMR(rfbE,raV,T2,muS,orbitType,Nrev,Ncase,LambOptions);

%orbit transfer 1 propagation
tspan_t1 = linspace( 0, T1, 1000 );
% manca viT1 velocità iniziale del primo transfer arc da salvare nella struttura results
[time_t1, state_t1] = ode113( @(t,s) tbp_ode(t, s, muS), tspan_t1, [rdJ, ViT1], options);

%orbit transfer 2 propagation
tspan_t2 = linspace( 0, T2, 1000 );
% manca viT2 velocità iniziale del secondo transfer arc da salvare nella struttura results
[time_t2, state_t2] = ode113( @(t,s) tbp_ode(t, s, muS), tspan_t2, [rfbE, ViT2], options);

figure(1)  
plot3( stateJ(:,1), stateJ(:,2), stateJ(:,3), 'red', 'linewidth', 2)        %Jupiter orbit
hold on
plot3( state_t1(:,1), state_t1(:,2), state_t1(:,3), 'green--', 'linewidth', 1) %transfer arc1
plot3( stateE(:,1), stateE(:,2), stateE(:,3), 'blue', 'linewidth', 2)       %Earth orbit
plot3( state_t2(:,1), state_t2(:,2), state_t2(:,3), 'magenta--', 'linewidth', 1) %transfer arc2 
plot3( stateV(:,1), stateV(:,2), stateV(:,3), 'yellow', 'linewidth', 2)       %Venus orbit


opts.units= 's';
planet3D('Sun',opts, [0 0 0]);  % add Planet

opts.units= 'j';
planet3D( 'Jupiter',opts, rdJ);

opts.units= 'p';
planet3D('Earth', opts, rfbE);
planet3D( 'Venus',opts, raV);

legend('Jupiter orbit','Transfer arc 1', 'Earth orbit','Transfer arc 2','Venus orbit', 'Sun', 'Jupiter at launch','Earth at flyby', 'Venus at arrival');

xlabel("r_x[km]");
ylabel("r_y[km]");
zlabel("r_z[km]");
title('Trajectory');
grid on
axis equal

%% flyby 

tspanH = linspace( 0, 1000, 1000 );
[timeH1, stateH1] = ode113( @(t,s) tbp_ode(t, s, muE), tspanH, [r_p(a,b,c,:),vp_p(a,b,c,:)], options); %outoging hyperbola
[timeH2, stateH2] = ode113( @(t,s) tbp_ode(t, s, muE), -tspanH, [r_p(a,b,c,:),vp_m(a,b,c,:)], options);%incoming hyperbola
r_H1= (stateH1(:,1:3));                            
r_H2= (stateH1(:,1:3));   

G= astroConstants(1);
ME= astroConstants(13)/G;
MS= astroConstants(4)/G;

r_soi= astroConstants(2)*(ME/MS)^(2/5);
flag1=1;
flag2=1;

for kk=1:length(tspanH)
    if norm(r_H1(kk,:)) > r_soi && flag1==1
        time1= timeH1(kk);
        flag1=0;
    end
    if norm(r_H2(kk,:)) > r_soi && flag2==1
        time2= -timeH2(kk);
        flag2=0;
    end
end

fb_time= time1+time2 % è un risutato importante va portato fuori come globale

[rp, r_p,Dv_p, vp_p, vp_m] = flyby(VfT1, ViT2, rfbE, vfbE);

tspanH = linspace( 0, 30, 1000 );
[timeH1, stateH1] = ode113( @(t,s) tbp_ode(t, s, muE), tspanH, [r_p,vp_p], options); %outoging hyperbola
[timeH2, stateH2] = ode113( @(t,s) tbp_ode(t, s, muE), -tspanH, [r_p,vp_m], options);%incoming hyperbola

figure(2)   %flyby in Earth SOI
p1=plot3( stateH1(:,1), stateH1(:,2), stateH1(:,3), 'red', 'linewidth', 2); %outgoing      
hold on
p2=plot3( stateH2(:,1), stateH2(:,2), stateH2(:,3), 'blue', 'linewidth', 2);  %incoming
opts.units= 'km';
planet3D('Earth',opts, [0 0 0]);  % add Planet
grid on
axis equal

%add arrow in direction of planet velocity
p3=arrow3d([0 vfbE(1)*10^3],[0 vfbE(2)*10^3],[0 vfbE(3)*10^3],9/10,0.2*10^3,0.8*10^3,'black');  
%add position of perigee
p4=scatter3(stateH1(1,1), stateH1(1,2),stateH1(1,3), 'blue' ,'filled', 'LineWidth',10);

title('Flyby in Earth-centred frame parallel to HECI')
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
legend([p2 p1 p3 p4],'incoming hyperbola', 'outgoing hyperbola', 'planet velocity', 'pericenter' )

%% preiminary analysis on single arc 1

in_dep_date= date2mjd2000([2025 08 01 00 00 00]);  %earliest launch date
fin_dep_date= date2mjd2000([2060 08 01 00 00 00]);  %latest launch date

in_fb_date= date2mjd2000([ 2025 08 01 00 00 00]);
fin_fb_date= date2mjd2000([ 2061 08 01 00 00 00]);

n=1000;

dep_window= linspace(in_dep_date, fin_dep_date,n);
fb_window= linspace(in_fb_date, fin_fb_date,n);

Dvmin=astroConstants(5);                                %initialize Vmin with a high value
r=NaN;
c=NaN;                                                  %initialize row and column of the min
muS= astroConstants(4);

piano= (1:1:5);                                    %plane to cut contour

options = odeset ('RelTol', 1e-3, 'AbsTol', 1e-3);

%initialize every matrix before the cycle to avoid warnings
r1= ones(3,n);
v1= ones(3,n);
r2= ones(3,n);
v2= ones(3,n);
T11= ones(n,n);
aT1= ones(n,n);
Vi1= ones(n,n,3);
Vf1= ones(n,n,3);
Dvtot1= ones(n,n);
x= ones(1,n);
y= ones(1,n);

for i=1:n
    kep1= uplanet(dep_window(i), 5);
    [r1(:,i), v1(:,i)]= kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),muS);

    for j=1:n
        kep2= uplanet(fb_window(j), 3);
        [r2(:,j), v2(:,j)]= kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muS);
        
        T11(i,j)= (fb_window(j)-dep_window(i))*24*60*60;    %ToF converted from Julian day to sec
        
        %Lambert prob relative to i-departure date and j-arrival date
        orbitType=0;                %direct orbit
        Nrev=0;                     %zero revolution case
        Ncase=0;
        LambOptions=1;
        [aT1(i,j),pT,E_T,ERROR,Vi1(i,j,:),Vf1(i,j,:),TPAR,THETA] = lambertMR(r1(:,i),r2(:,j),T11(i,j),muS,orbitType,Nrev,Ncase,LambOptions);
        Dv1x= Vi1(i,j,1)-v1(1,i);
        Dv1y= Vi1(i,j,2)-v1(2,i);
        Dv1z= Vi1(i,j,3)-v1(3,i);
        Dv1= norm([Dv1x Dv1y Dv1z]);    %delta V of manoeuvre 1 (from orbit 1 to transfer
        Dvtot(i,j)= Dv1;    
        %research for min values of Dvtot
        if Dvtot(i,j)<Dvmin 
            Dvmin= Dvtot(i,j);
            r=i;
            c=j;
        end
    end
    x(i)=datenum(mjd20002date(dep_window(i)));  %conversion julian day to regular date to matlab date form
    y(i)=datenum(mjd20002date(fb_window(i)));  %conversion julian day to regular date to matlab date form
end 
   launch_day= mjd20002date(dep_window(r));
   flyby_day= mjd20002date(fb_window(c));


%% contour pork chop plot
figure(1)
hold on
[C, h] = contour(x, y, Dvtot', floor(Dvmin) + (piano));
plot(x(r),y(c), '*', 'LineWidth',20);
caxis([Dvmin Dvmin + piano(length(piano))]);                                     % Cambiare il piano di sezione in un intorno del minimo [in particolare, sopra il MIN]
Color= colorbar;                                              % Aggiungere la barra dei colori a destra
clabel(C,h, floor(Dvmin) + (piano));                          % Aggiungere le etichette sulle curve di livello
datetick( 'x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
datetick( 'y', 'yyyy mmm dd', 'keeplimits');
xtickangle(45);
xlabel('Jupiter departure date');
ylabel('Earth flyby date');
set(Color,'YtickLabel');
Color.Label.String= 'Delta V1 [km/s]';
grid on

title('Pork Chop Plot');

[C1, h1] = contour(x, y, (T11/(24*3600))','Black');
clabel(C1,h1);                                               % Aggiungere le etichette sulle curve di livello

%% plot trajectory

TJ = 2*pi*sqrt( (kep1(1))^3/muS );              %period of Jupiter around the Sun
tspanJ = linspace( 0, TJ, 1000 );
[timeJ, stateJ] = ode113( @(t,s) tbp_ode(t, s, muS), tspanJ, [r1(:,r), v1(:,r)], options);

TE = 2*pi*sqrt( (kep2(1))^3/muS );              %period of Earth around the Sun
tspanE = linspace( 0, TE, 1000 );
[timeE, stateE] = ode113( @(t,s) tbp_ode(t, s, muS), tspanE, [r2(:,c), v2(:,c)], options);

%orbit transfer propagation
tspan_t1 = linspace( 0, T11(r,c), 1000 );
[time_t1, state_t1] = ode113( @(t,s) tbp_ode(t, s, muS), tspan_t1, [r1(:,r), [Vi1(r,c,1) Vi1(r,c,2) Vi1(r,c,3)]'], options);

figure(2)          
plot3( stateJ(:,1), stateJ(:,2), stateJ(:,3), 'blue', 'linewidth', 2)       %Jupiter orbit
hold on
plot3( state_t1(:,1), state_t1(:,2), state_t1(:,3), 'green--', 'linewidth', 2) %transfer arc
plot3( stateE(:,1), stateE(:,2), stateE(:,3), 'red', 'linewidth', 2)        %Earth orbit

opts.units= 's';
planet3D('Sun',opts, [0 0 0]);  % add Planet

opts.units= 'p';
planet3D('Jupiter', opts, r1(:,r));
planet3D( 'Earth',opts, r2(:,c));

legend('Jupiter orbit','Transfer arc 1', 'Earth orbit', 'Sun', 'Jupiter at launch', 'Earth at flyby');

xlabel("r_x[km]");
ylabel("r_y[km]");
zlabel("r_z[km]");
title('Trajectory');
grid on
axis equal

%% %% preiminary analysis on single arc 2

in_arr_date= date2mjd2000([2027 08 01 00 00 00]);
%given
fin_arr_date= date2mjd2000([2065 08 01 00 00 00]);  

%tof ottimale inferiore a 1000 giorni

n=1000;      %has to be picked reasoning
arr_window= linspace(in_arr_date, fin_arr_date,n);


Dvmin=astroConstants(5);                                %initialize Vmin with a high value
r=NaN;
c=NaN;                                                  %initialize row and column of the min
muS= astroConstants(4);

piano= (1:1:10);                                    %plane to cut contour

options = odeset ('RelTol', 1e-3, 'AbsTol', 1e-3);

%initialize every matrix before the cycle to avoid warnings
ra3= ones(3,n);
va3= ones(3,n);
r2= ones(3,n);
v2= ones(3,n);
T22= ones(n,n);
aT2= ones(n,n);
Vi2= ones(n,n,3);
Vf2= ones(n,n,3);
Dvtot= ones(n,n);
x= ones(1,n);
y= ones(1,n);

for j=1:n
    kep2= uplanet(fb_window(j), 3);
    [r2(:,j), v2(:,j)]= kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muS);

    for k=1:n
        kep3= uplanet(arr_window(k), 2);
        [r3(:,k), v3(:,k)]= kep2car(kep3(1),kep3(2),kep3(3),kep3(4),kep3(5),kep3(6),muS);
        
        T22(j,k)= (arr_window(k)-fb_window(j))*24*60*60;    %ToF converted from Julian day to sec
        
        %Lambert prob relative to j-departure date and k-arrival date
        orbitType=0;                %direct orbit
        Nrev=0;                     %zero revolution case
        Ncase=0;
        LambOptions=1;
        [aT2(j,k),pT,E_T,ERROR,Vi2(j,k,:),Vf2(j,k,:),TPAR,THETA] = lambertMR(r2(:,j),r3(:,k),T22(j,k),muS,orbitType,Nrev,Ncase,LambOptions);
        Dv2x= v3(1,k)-Vf2(j,k,1);
        Dv2y= v3(2,k)-Vf2(j,k,2);
        Dv2z= v3(3,k)-Vf2(j,k,3);
        Dv2= norm([Dv2x Dv2y Dv2z]);    %delta V of manoeuvre 2 (from transfer 2 to Venus orbit)
        
        Dvtot(j,k)= Dv2;    
        %research for min values of Dvtot
        if Dvtot(j,k)<Dvmin 
            Dvmin= Dvtot(j,k);
            r=j;
            c=k;
        end
    end
    x(j)=datenum(mjd20002date(fb_window(j)));  %conversion julian day to regular date to matlab date form
    y(j)=datenum(mjd20002date(arr_window(j)));  %conversion julian day to regular date to matlab date form
end 
   flyby_day= mjd20002date(fb_window(r));
   arrival_day= mjd20002date(arr_window(c));


%% contour pork chop plot
figure(1)
hold on
[C, h] = contour(x, y, Dvtot', floor(Dvmin) + (piano));
plot(x(r),y(c), '*', 'LineWidth',20);
caxis([Dvmin Dvmin + piano(length(piano))]);                                     % Cambiare il piano di sezione in un intorno del minimo [in particolare, sopra il MIN]
Color= colorbar;                                              % Aggiungere la barra dei colori a destra
clabel(C,h, floor(Dvmin) + (piano));                          % Aggiungere le etichette sulle curve di livello
datetick( 'x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
datetick( 'y', 'yyyy mmm dd', 'keeplimits');
xtickangle(45);
xlabel('Earth flyby date');
ylabel('Venus arrival date');
set(Color,'YtickLabel');
Color.Label.String= 'Delta V2 [km/s]';
grid on

title('Pork Chop Plot');

[C1, h1] = contour(x, y, (T22/(24*3600))','Black');
clabel(C1,h1);                                               % Aggiungere le etichette sulle curve di livello

%% plot trajectory

TE = 2*pi*sqrt( (kep2(1))^3/muS );              %period of Earth around the Sun
tspanE = linspace( 0, TE, 1000 );
[timeE, stateE] = ode113( @(t,s) tbp_ode(t, s, muS), tspanE, [r2(:,r), v2(:,r)], options);

TV = 2*pi*sqrt( (kep3(1))^3/muS );              %period of Earth around the Sun
tspanV = linspace( 0, TV, 1000 );
[timeV, stateV] = ode113( @(t,s) tbp_ode(t, s, muS), tspanV, [r3(:,c), v3(:,c)], options);

%orbit transfer propagation
tspan_t2 = linspace( 0, T22(r,c), 1000 );
[time_t2, state_t2] = ode113( @(t,s) tbp_ode(t, s, muS), tspan_t2, [r2(:,r), [Vi2(r,c,1) Vi2(r,c,2) Vi2(r,c,3)]'], options);

figure(2)          
plot3( stateE(:,1), stateE(:,2), stateE(:,3), 'blue', 'linewidth', 2)       %Earth orbit
hold on
plot3( state_t2(:,1), state_t2(:,2), state_t2(:,3), 'green--', 'linewidth', 2) %transfer arc 2
plot3( stateV(:,1), stateV(:,2), stateV(:,3), 'red', 'linewidth', 2)        %Earth orbit

opts.units= 's';
planet3D('Sun',opts, [0 0 0]);  % add Planet

opts.units= 'p';
planet3D('Earth', opts, r2(:,r));
planet3D( 'Venus',opts, r3(:,c));

legend('Earth orbit','Transfer arc 2', 'Venus orbit', 'Sun', 'Earth at flyby', 'Venus at arrival');

xlabel("r_x[km]");
ylabel("r_y[km]");
zlabel("r_z[km]");
title('Trajectory');
grid on
axis equal
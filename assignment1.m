%% project assignment 1

% group id 2163
% departure planet Jupiter
% flyby Earth
% arrival planet Venus
% earliest departure 2025/08/01
% latest arrival     2065/08/01

clearvars
close all 
clc

muE= astroConstants(13);
muS= astroConstants(4);
muJ= astroConstants(15);
muV= astroConstants(12);

RE= astroConstants(23);              %Earth radius
h_atm= 100;                          %height Earth atmosphere in km

%given
in_dep_date= date2mjd2000([2035 11 01 00 00 00]);  %earliest launch date
fin_arr_date= date2mjd2000([2040 08 01 00 00 00]);  %latest arrival date
% have to be picked reasoning how and why
%picked randomly right now
in_fb_date= date2mjd2000([ 2038 05 01 00 00 00]);
fin_fb_date= date2mjd2000([ 2038 07 01 00 00 00]);
fin_dep_date= date2mjd2000([2036 01 01 00 00 00]);  %latest launch date
in_arr_date= date2mjd2000([2037 01 01 00 00 00]);   %earliest arrival date

n=30;      %has to be picked reasoning
%can and should be different for dep_window, fb_window and arr_window

tof1= linspace(0,2000*24*60*60, n);
tof2= linspace(0,1000*24*60*60, n);

dep_window= linspace(in_dep_date, fin_dep_date,n);
fb_window= linspace(in_fb_date, fin_fb_date,n);
arr_window= linspace(in_arr_date, fin_arr_date,n);

Dvmin=astroConstants(5);                                %initialize Vmin with a high value

%initialize variables 
rdJ= ones(3,n);
vdJ= ones(3,n);
rfbE= ones(3,n);
vfbE= ones(3,n);
raV= ones(3,n);
vaV= ones(3,n);

T1= ones(n,n);
aT1= ones(n,n);
ViT1= ones(n,n,3);
VfT1= ones(n,n,3);
DvT1= ones(n,n);

T2= ones(n,n);
aT2= ones(n,n);
ViT2= ones(n,n,3);
VfT2= ones(n,n,3);
DvT2= ones(n,n);

DviT1=ones(n,n);
DvfT2=ones(n,n);
Dvtot= ones(n,n,n);
%Tof=ones(n,n,n);

for i= 1:n
    kepJ= uplanet(dep_window(i), 5);
    [rdJ(:,i), vdJ(:,i)]= kep2car(kepJ(1),kepJ(2),kepJ(3),kepJ(4),kepJ(5),kepJ(6),muS);

    for j=1:n
        kepE= uplanet(fb_window(j), 3);
        [rfbE(:,j), vfbE(:,j)]= kep2car(kepE(1),kepE(2),kepE(3),kepE(4),kepE(5),kepE(6),muS);
        T1(i,j)= (fb_window(j)-dep_window(i))*24*60*60;     %time of flight Jupiter-Earth

        %solve Lambert arc 1
        orbitType=0;                %direct orbit
        Nrev=0;                     %zero revolution case
        Ncase=0;
        LambOptions=1;
        [aT1(i,j),pT,E_T,ERROR,ViT1(i,j,:),VfT1(i,j,:),TPAR,THETA] = lambertMR(rdJ(:,i),rfbE(:,j),T1(i,j),muS,orbitType,Nrev,Ncase,LambOptions);

        DviT1x= ViT1(i,j,1)-vdJ(1,i);
        DviT1y= ViT1(i,j,2)-vdJ(2,i);
        DviT1z= ViT1(i,j,3)-vdJ(3,i);
        DviT1(i,j)= norm([DviT1x DviT1y DviT1z]);    %delta V of manoeuvre 1 (from Jupiter's orbit to transfer1)
       

        DvT1(i,j)= DviT1(i,j); 

        for k=1:n
            kepV= uplanet(arr_window(k), 2);
            [raV(:,k), vaV(:,k)]= kep2car(kepV(1),kepV(2),kepV(3),kepV(4),kepV(5),kepV(6),muS);
    
            T2(j,k)= (arr_window(k)-fb_window(j))*24*60*60;        %time of flight Earth-Venus
            %ToF(i,j,k)= arr_window(k)-dep_window(i);    %total time of flight Juiter-Venus

            %solve the relative Lambert arc
            orbitType=0;                %direct orbit
            Nrev=0;                     %zero revolution case
            Ncase=0;
            LambOptions=1;
            [aT2(j,k),pT,E_T,ERROR,ViT2(j,k,:),VfT2(j,k,:),TPAR,THETA] = lambertMR(rfbE(:,j),raV(:,k),T2(j,k),muS,orbitType,Nrev,Ncase,LambOptions);
       
            DvfT2x= vaV(1,k)-VfT2(j,k,1);
            DvfT2y= vaV(2,k)-VfT2(j,k,2);
            DvfT2z= vaV(3,k)-VfT2(j,k,3);
            DvfT2(j,k)= norm([DvfT2x DvfT2y DvfT2z]);

            DvT2(j,k)= DvfT2(j,k);

            % Earth flyby
            [rp, r_p(i,j,k,:),Dv_p(i,j,k), vp_p(i,j,k,:), vp_m(i,j,k,:)] = flyby(VfT1(i,j,:), ViT2(j,k,:), rfbE(:,j), vfbE(:,j));
            
            Dvtot(i,j,k)= DvT1(i,j) + DvT2(j,k)+ Dv_p(i,j,k) ;

            %research for min values of Dvtot
            if Dvtot(i,j,k)<Dvmin && rp> RE+h_atm
              Dvmin= Dvtot(i,j,k);
              %position of the minimum
              a=i;  
              b=j;
              c=k;
              h_fb= rp-RE;    %alttude of the flyby
            end 
        end
    end
end

%print useful results 
DeltaV_totale= Dvtot(a,b,c)
DeltaV_arc1= DvT1(a,b)
DeltaV_arc2= DvT2(b,c)
DeltaV_flyby_powered= Dv_p(a,b,c)
DeltaV_flyby_natural= ViT2(b,c,:)-VfT1(a,b,:)

launch_day= mjd20002date(dep_window(a))
flyby_day= mjd20002date(fb_window(b))
arrival_day= mjd20002date(arr_window(c))

%% plot heliocnetric trajectory

options = odeset ('RelTol', 1e-3, 'AbsTol', 1e-3);

TJ = 2*pi*sqrt( (kepJ(1))^3/muS );              %period of Earth around the Sun
tspanJ = linspace( 0, TJ, 1000 );
[timeJ, stateJ] = ode113( @(t,s) tbp_ode(t, s, muS), tspanJ, [rdJ(:,a), vdJ(:,a)], options);

TE = 2*pi*sqrt( (kepE(1))^3/muS );              %period of Earth around the Sun
tspanE = linspace( 0, TE, 1000 );
[timeE, stateE] = ode113( @(t,s) tbp_ode(t, s, muS), tspanE, [rfbE(:,b), vfbE(:,b)], options);

TV = 2*pi*sqrt( (kepV(1))^3/muS );              %period of Earth around the Sun
tspanV = linspace( 0, TV, 1000 );
[timeV, stateV] = ode113( @(t,s) tbp_ode(t, s, muS), tspanV, [raV(:,c), vaV(:,c)], options);

%orbit transfer 1 propagation
tspan_t1 = linspace( 0, T1(a,b), 1000 );
[time_t1, state_t1] = ode113( @(t,s) tbp_ode(t, s, muS), tspan_t1, [rdJ(:,a), [ViT1(a,b,1) ViT1(a,b,2) ViT1(a,b,3)]'], options);

%orbit transfer 2 propagation
tspan_t2 = linspace( 0, T2(b,c), 1000 );
[time_t2, state_t2] = ode113( @(t,s) tbp_ode(t, s, muS), tspan_t2, [rfbE(:,b), [ViT2(b,c,1) ViT2(b,c,2) ViT2(b,c,3)]'], options);



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
planet3D( 'Jupiter',opts, rdJ(:,a));
%planet3D( 'Jupiter',opts, rdJ(:,b));
%planet3D( 'Jupiter',opts, rdJ(:,c));

opts.units= 'p';
planet3D('Earth', opts, rfbE(:,a));
planet3D('Earth', opts, rfbE(:,b));
%planet3D('Earth', opts, rfbE(:,c));

planet3D( 'Venus',opts, raV(:,a));
planet3D( 'Venus',opts, raV(:,b));
planet3D( 'Venus',opts, raV(:,c));


legend('Jupiter orbit','Transfer arc 1', 'Earth orbit','Transfer arc 2','Venus orbit', 'Sun', 'Jupiter at launch','Earth at launch','Earth at flyby','Venus at launch', 'Venus at flyby', 'Venus at arrival');

xlabel("r_x[km]");
ylabel("r_y[km]");
zlabel("r_z[km]");
title('Trajectory');
grid on
axis equal

%% plot flyby
 
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

fb_time= time1+time2
%%
tspanH = linspace( 0, 30, 1000 );
[timeH1, stateH1] = ode113( @(t,s) tbp_ode(t, s, muE), tspanH, [r_p(a,b,c,:),vp_p(a,b,c,:)], options); %outoging hyperbola
[timeH2, stateH2] = ode113( @(t,s) tbp_ode(t, s, muE), -tspanH, [r_p(a,b,c,:),vp_m(a,b,c,:)], options);%incoming hyperbola

figure(2)   %flyby in Earth SOI
p1=plot3( stateH1(:,1), stateH1(:,2), stateH1(:,3), 'red', 'linewidth', 2); %outgoing      
hold on
p2=plot3( stateH2(:,1), stateH2(:,2), stateH2(:,3), 'blue', 'linewidth', 2);  %incoming
opts.units= 'km';
planet3D('Earth',opts, [0 0 0]);  % add Planet
grid on
axis equal

%add arrow in direction of planet velocity
p3=arrow3d([0 vfbE(1,b)*10^3],[0 vfbE(2,b)*10^3],[0 vfbE(3,b)*10^3],9/10,0.2*10^3,0.8*10^3,'black');  
%add position of perigee
%p4=scatter3(r_p(a,b,c,1), r_p(a,b,c,2), r_p(a,b,c,3), 'blue' ,'filled', 'LineWidth',10);
p4=scatter3(stateH1(1,1), stateH1(1,2),stateH1(1,3), 'blue' ,'filled', 'LineWidth',10);

title('Flyby in Earth-centred frame parallel to HECI')
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')
legend([p2 p1 p3 p4],'incoming hyperbola', 'outgoing hyperbola', 'planet velocity', 'pericenter' )
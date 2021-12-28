%% assignement 2

run config.m

% satellite period
T = 2*pi/sqrt(muE/a0^3); % [s]

% cartesian initial conditions
[r0 v0] = kep2car(a0, e0, deg2rad(i0), deg2rad(Ome0), deg2rad(ome0), deg2rad(f0), muE);

% time interval of integration 
tspan = linspace(MJD2000*24*3600, MJD2000*24*3600 + T, 1000);

% integration
[time state] = ode113( @(t,s) tbp_ode(t,s,muE), tspan, [r0; v0], options_sol);

% plotting the nominal orbit
earth_sphere
hold on 
plot3( state(:,1), state(:,2), state(:,3), 'r', 'LineWidth', 2 )
xlabel('r_x [km]')
ylabel('r_y [km]')
zlabel('r_z [km]')
title ('Nominal orbit in an unperturbed condition') 
axis equal  
grid on

%% Groundtracks

% calculation of a for repeating groundtrack unperturbed
a_rep_unp = repeatingGroundTracks(k, m, muE, omeE);

Kep_elements_rep_unp = [a_rep_unp, e0, deg2rad(i0), deg2rad(Ome0), deg2rad(ome0), deg2rad(f0)];

x = input('Enter an interval of time: ');

switch x
    case '1 period'
        tspan = linspace(MJD2000*24*3600, MJD2000*24*3600 + T, 100000);
        b = 'one period';

    case '1 day'
        tspan = linspace(MJD2000*24*3600, (MJD2000+1)*3600*24, 100000);
        b = 'one day';

    case '10 days'
        tspan = linspace(MJD2000*24*3600, (MJD2000+10)*3600*24,100000);
        b = '10 days';
end

[alpha, delta, lon, lat] = groundTrack (Kep_elements, deg2rad(lon0), tspan, omeE, muE, t0);
[alpha_rep_unp, delta_rep_unp, lon_rep_unp, lat_rep_unp] = groundTrack (Kep_elements_rep_unp, deg2rad(lon0), tspan, omeE, muE, t0);
[alpha_per, delta_per, lon_per, lat_per] = groundTrack_perturbed(Kep_elements, deg2rad(lon0), tspan, omeE, muE, Re, J2, cr, A2m, t0);
[alpha_per_rep, delta_per_rep, lon_per_rep, lat_per_rep] = groundTrack_perturbed(Kep_elements_rep_unp, deg2rad(lon0), tspan, omeE, muE, Re, J2, cr, A2m, t0);


figure()

% to plot the earth as background
image_file = 'earth_2D.jpg'; % saved in the same folder
cdata = imread(image_file);
image([-180 180],[90 -90], cdata);
ax = gca;
ax.YDir = 'normal'; % to rotate the earth draw

% fig1 = figure('Name','Ground Track 1 orbit');
% set(fig1, 'Units', 'Normalized', 'OuterPosition', [.15 .25 .7 .7]);
% hold on;
% axis equal;
% set(gca,'XTick',[-180:30:180],'XTickMode','manual');
% set(gca,'YTick',[-90:30:90],'YTickMode','manual');
% xlim([-180,180]); ylim([-90,90]);
%         
% image_file = 'earth_2D.jpg';
% cdata      = flip(imread(image_file));
% imagesc([-180,180],[-90, 90],cdata);

hold on 

% plot of unperturbed case
plot( wrapTo180(rad2deg(lon)), rad2deg(lat),'.', 'Color','g');
plot( wrapTo180(rad2deg(lon(1))), rad2deg(lat(1)),'*', 'LineWidth', 10, 'Color', 'g');
plot( wrapTo180(rad2deg(lon(length(lon)))), rad2deg(lat(length(lat))),'^', 'LineWidth', 6, 'Color', 'g');


% plot of repeating unperturbed case
plot( wrapTo180(rad2deg(lon_rep_unp)), rad2deg(lat_rep_unp),'.', 'Color','r');
plot( wrapTo180(rad2deg(lon_rep_unp(1))), rad2deg(lat_rep_unp(1)),'*', 'LineWidth', 10, 'Color', 'r');
plot( wrapTo180(rad2deg(lon_rep_unp(length(lon_rep_unp)))), rad2deg(lat_rep_unp(length(lat_rep_unp))),'^', 'LineWidth', 6, 'Color', 'r');

% plot of perturbed case (J2 + SRP)
plot( wrapTo180(rad2deg(lon_per)), rad2deg(lat_per),'.', 'LineWidth', 7, 'Color','b');
plot( wrapTo180(rad2deg(lon_per(1))), rad2deg(lat_per(1)),'*', 'LineWidth', 10, 'Color', 'b');
plot( wrapTo180(rad2deg(lon_per(length(lon_per)))), rad2deg(lat_per(length(lat_per))),'^', 'LineWidth', 6, 'Color', 'b');

legend('Unperturbed', 'start unperturbed', 'end unperturbed', 'Repeating unperturbed', 'start repeating unperturbed', 'end repeating unperturbed', 'Perturbed', 'start perturbed', 'end perturbed');

axis ([-180 180 -90 90]);
axis equal;
grid on;
xlabel('lon [deg]');
ylabel('lat [deg]');
% title("Groundtracks case "b" ");

figure()

% to plot the earth as background
image_file = 'earth_2D.jpg'; % saved in the same folder
cdata = imread(image_file);
image([-180 180],[90 -90], cdata);
ax = gca;
ax.YDir = 'normal'; % to rotate the earth draw

hold on

% plot of repeating unperturbed case
plot( wrapTo180(rad2deg(lon_rep_unp)), rad2deg(lat_rep_unp),'.', 'Color','r');
plot( wrapTo180(rad2deg(lon_rep_unp(1))), rad2deg(lat_rep_unp(1)),'*', 'LineWidth', 10, 'Color', 'r');
plot( wrapTo180(rad2deg(lon_rep_unp(length(lon_rep_unp)))), rad2deg(lat_rep_unp(length(lat_rep_unp))),'^', 'LineWidth', 6, 'Color', 'r');

% plot of perturbed case with the modified semimajor axis for point b
plot( wrapTo180(rad2deg(lon_per_rep)), rad2deg(lat_per_rep),'.', 'LineWidth', 7, 'Color','y');
plot( wrapTo180(rad2deg(lon_per(1))), rad2deg(lat_per_rep(1)),'*', 'LineWidth', 10, 'Color', 'y');
plot( wrapTo180(rad2deg(lon_per_rep(length(lon_per_rep)))), rad2deg(lat_per_rep(length(lat_per_rep))),'^', 'LineWidth', 6, 'Color', 'y');

legend('Repeating unperturbed', 'start repeating unperturbed', 'end repeating unperturbed', 'Repeating perturbed', 'start repeating perturbed', 'end repeating perturbed');

axis ([-180 180 -90 90]);
axis equal;
grid on;
xlabel('lon [deg]');
ylabel('lat [deg]');
% title("Groundtracks case "b" ");


%% Introduction of the assigned perturbations J2 + SRP

% initial keplerian elements
kep0 = [a0 e0 deg2rad(i0) deg2rad(Ome0) deg2rad(ome0) deg2rad(f0)];

% time interval of integration 
tspan = linspace(MJD2000*24*3600, MJD2000*24*3600 + 3*365*3600*24, 40000);

% integration solver options
options_sol = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14); 

% propagation of Cartesian elements
[ time_car, state ] = ode113( @(t,s) ode_perturbed (t, s, muE,  Re, J2, cr, A2m), tspan, [r0; v0], options_sol ); 

% propagation of Keplerian elements
[ time_Kep, Kep ] = ode113( @(t,s) gauss_propagation(t, s, muE, Re, J2, cr, A2m), tspan, kep0', options_sol ); 

% Keplerian elements to cartesian elements
for j = 1:length(tspan)
    [state_Kep(j,1:3), state_Kep(j, 4:6)] = kep2car(Kep(j,1), Kep(j,2), Kep(j,3), Kep(j,4), Kep(j,5), Kep(j,6), muE);

end

figure()
% propagation of the orbit in Cartesian elements
plot3(state(:,1), state(:,2), state(:,3));
hold on 
earth_sphere
xlabel('r_x [km]')
ylabel('r_y [km]')
zlabel('r_z [km]')
axis equal
grid on
title("Propagation of the perturbed orbit in Cartesian elements");

figure()
% propagation of the orbit in Cartesian elements from Keplerian elements
% converted
plot3(state_Kep(:,1), state_Kep(:,2), state_Kep(:,3));
hold on 
earth_sphere
xlabel('r_x [km]')
ylabel('r_y [km]')
zlabel('r_z [km]')
axis equal
grid on
title("Propagation of the perturbed orbit in Cartesian elements from Keplerian elements");

%% propagation of Keplerian elements

% Keplerian elements to cartesian elements
for k = 1:length(tspan)
    [ Kep_state(k,1),  Kep_state(k,2),  Kep_state(k,3),  Kep_state(k,4),  Kep_state(k,5),  Kep_state(k,6) ] = ...
        car2kep(state(k,1:3), state(k,4:6), muE);
end

% to make RAAD, argument of perigee and true anomaly continuous
Kep_state(:,6) = unwrap(Kep_state(:,6));
Kep_state(:,4) = unwrap(Kep_state(:,4));
Kep_state(:,5) = unwrap(Kep_state(:,5));

y = input('Enter a or e: ');

switch y
    case 'a'
        name = 'a';
        par = 1;
        norm = kep0(1);
    
    case 'e'
        name = 'e';
        par = 2;
        norm = 1;
        
end

figure()
plot ((tspan-MJD2000*24*3600)/T, abs(Kep(:,par) - Kep_state(:,par))/norm);
grid on
% title ("relative error of 'name' between gauss and cartesian propagation")

figure()
plot ((tspan-MJD2000*24*3600)/T, Kep(:,par));
hold on
plot ((tspan-MJD2000*24*3600)/T, Kep_state(:,par));
grid on
legend ('Gauss', 'Cartesian')
% title ("propagation of 'name' ")


z = input('Enter i, Om, om or f: ');

switch z
     case 'i'
        name = 'i';
        par =3;
        norm = 2*pi;

    case 'Om'
        name = 'Om';
        par = 4;
        norm = 2*pi;

    case 'om'
        name = 'om';
        par = 5;
        norm = 2*pi;

    case 'f'
        name = 'f';
        par = 6;
        norm = Kep_state(:,6);

end    

figure()
plot ((tspan-MJD2000*24*3600)/T, rad2deg(abs(Kep(:,par) - Kep_state(:,par)))./abs(norm));
grid on
titlestr = strcat("Relative error of", name, "between gauss and cartesian propagation")
title (titlestr)

figure()
plot ((tspan-MJD2000*24*3600)/T, rad2deg(Kep(:,par)));
hold on
plot ((tspan-MJD2000*24*3600)/T, rad2deg(Kep_state(:,par)));
grid on
legend ('Gauss', 'Cartesian')
% title ("propagation of 'name' ")

%% real obervation
tspan_real = linspace(MJD2000*24*3600, MJD2000*24*3600 + 3*365*3600*24, length(A)/2);
tspan = linspace(MJD2000*24*3600, MJD2000*24*3600 + 3*365*3600*24, 40000);

figure()
% plot ((tspan-MJD2000*24*3600)/T, Kep(:,1));
% hold on
plot ((tspan-MJD2000*24*3600)/T, Kep(:,1), (tspan_real-MJD2000*24*3600)/T, Kep_real(:,1));
grid on
legend ('Gauss', 'Real')

%% filtering Keplerian elements

x = input('Enter a or e to be filtered: ');

switch x
    case 'a'
        name = 'a';
        par = 1;
        % number of elements in the movemean
        k_par = T/3600/T*length(tspan);
    
    case 'e'
        name = 'e';
        par = 2;
        % number of elements in the movemean
        k_par = T/3600/T*length(tspan);
%         k_par2 = 800*T/3600/T*length(tspan);

end

Kep_filt = movmean(Kep(:,par), k_par);
% Kep_filt2 = movmean(Kep(:,par), k_par2);

figure()
plot ((tspan-MJD2000*24*3600)/T, Kep(:,par));
hold on
plot ((tspan-MJD2000*24*3600)/T, Kep_filt);
% plot ((tspan-MJD2000*24*3600)/T, Kep_filt2);
grid on
legend ('Original', 'Filtered')

x = input('Enter i, Om, om or f to be filtered: ');

switch x
    case 'i'
        name = 'i';
        par = 3;
        % number of elements in the movemean
        k_par = T/3600/T*length(tspan);
%       k_par2 = 670*T/3600/T*length(tspan);
    
    case 'Om'
        name = 'Om';
        par = 4;
        % number of elements in the movemean
        k_par = T/3600/T*length(tspan);

        case 'om'
        name = 'om';
        par = 5;
        % number of elements in the movemean
        k_par = T/3600/T*length(tspan);

        case 'f'
        name = 'f';
        par = 6;
        % number of elements in the movemean
        k_par = T/3600/T*length(tspan);
end

Kep_filt = movmean(Kep(:,par), k_par);
% Kep_filt1 = movmean(Kep(:,par), k_par1);
% Kep_filt2 = movmean(Kep(:,par), k_par2);

figure()
plot ((tspan-MJD2000*24*3600)/T, rad2deg(Kep(:,par)));
hold on
plot ((tspan-MJD2000*24*3600)/T, rad2deg(Kep_filt));
grid on
legend ('Original', 'Filtered')

% figure()
% plot ((tspan-MJD2000*24*3600)/T, Kep(:,par));
% hold on
% plot ((tspan-MJD2000*24*3600)/T, Kep_filt);
% plot ((tspan-MJD2000*24*3600)/T, Kep_filt2);
% grid on
% legend ('Original', 'Filtered', 'Filtered1')
% 
% 

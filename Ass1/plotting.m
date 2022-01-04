% PLOTTING FILE
% In this file all graphs are plotted

%-------------------------------------------------------------------------%

%% PLOT INTERPLANETARY TRAJECTORY

% Jupiter's position at selected departure date
kepJ = uplanet(results.dates(1), 5);
[rdJ, vdJ] = kep2car(kepJ(1),kepJ(2),kepJ(3),kepJ(4),kepJ(5),kepJ(6),astro.muS);
% Jupiter's position at selected flyby date
kepJ2 = uplanet(results.dates(2), 5);
[rfbJ, vfbJ2] = kep2car(kepJ2(1),kepJ2(2),kepJ2(3),kepJ2(4),kepJ2(5),kepJ2(6),astro.muS);
% Jupiter's position at selected arrival date
kepJ3 = uplanet(results.dates(3), 5);
[raJ, vaJ] = kep2car(kepJ3(1),kepJ3(2),kepJ3(3),kepJ3(4),kepJ3(5),kepJ3(6),astro.muS);

% Earth's position at selected departure date
kepE1 = uplanet(results.dates(1), 3);
[rdE, vdE] = kep2car(kepE1(1),kepE1(2),kepE1(3),kepE1(4),kepE1(5),kepE1(6),astro.muS);
% Earth's position at selected flyby date
kepE = uplanet(results.dates(2), 3);
[rfbE, vfbE] = kep2car(kepE(1),kepE(2),kepE(3),kepE(4),kepE(5),kepE(6),astro.muS);
% Earth's position at selected arrival date
kepE3 = uplanet(results.dates(3), 3);
[raE, vaE] = kep2car(kepE3(1),kepE3(2),kepE3(3),kepE3(4),kepE3(5),kepE3(6),astro.muS);

% Venus' position at selected departure date
kepV1 = uplanet(results.dates(1), 2);
[rdV, vdV] = kep2car(kepV1(1),kepV1(2),kepV1(3),kepV1(4),kepV1(5),kepV1(6),astro.muS);
% Venus' position at selected flyby date
kepV2 = uplanet(results.dates(2), 2);
[rfbV, vfbV] = kep2car(kepV2(1),kepV2(2),kepV2(3),kepV2(4),kepV2(5),kepV2(6),astro.muS);
% Venus' position at selected arrival date
kepV = uplanet(results.dates(3), 2);
[raV, vaV] = kep2car(kepV(1),kepV(2),kepV(3),kepV(4),kepV(5),kepV(6),astro.muS);

% Jupiter's orbit propagation
TJ = 2*pi*sqrt((kepJ(1))^3/astro.muS );              
tspanJ = linspace(0, TJ, temp.n);
[~, stateJ] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspanJ, [rdJ, vdJ], options.ode_options);
% Earth's orbit propagation
TE = 2*pi*sqrt((kepE(1))^3/astro.muS );             
tspanE = linspace(0, TE, temp.n);
[~, stateE] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspanE, [rfbE, vfbE], options.ode_options);
% Venus' orbit propagation
TV = 2*pi*sqrt((kepV(1))^3/astro.muS );              
tspanV = linspace(0, TV, temp.n);
[~, stateV] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspanV, [raV, vaV], options.ode_options);

% Lambert arc 1 and 2 propagation
% 1st Lambert arc delta time
T1 = (results.dates(2)-results.dates(1))*24*60*60;
% 1st Lambert arc delta time
T2 = (results.dates(3)-results.dates(2))*24*60*60;

% 1st Lambert arc propagation
tspan_t1 = linspace(0, T1, temp.n);
[~, state_t1] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspan_t1, [rdJ, results.ViT1'], options.ode_options);
% 2nd Lambert arc propagation
tspan_t2 = linspace(0, T2, temp.n);
[~, state_t2] = ode113( @(t,s) tbp_ode(t, s, astro.muS), tspan_t2, [rfbE, results.ViT2'], options.ode_options);

flyby = figure();
flyby.WindowState = 'maximized';
view(0,90)
plot3( stateJ(:,1), stateJ(:,2), stateJ(:,3), 'red:', 'linewidth', 1.5)  
hold on
plot3( state_t1(:,1), state_t1(:,2), state_t1(:,3), 'green', 'linewidth', 1)
plot3( stateE(:,1), stateE(:,2), stateE(:,3), 'blue:', 'linewidth', 1.5)     
plot3( state_t2(:,1), state_t2(:,2), state_t2(:,3), 'magenta', 'linewidth', 1) 
plot3( stateV(:,1), stateV(:,2), stateV(:,3), 'yellow:', 'linewidth', 1.5)    

% add Planets at corresponding positions
opts.units = 's';
planet3D('Sun', opts, [0 0 0]);

opts.units = 'j';
planet3D( 'Jupiter', opts, rdJ);
planet3D( 'Jupiter', opts, rfbJ);
planet3D( 'Jupiter', opts, raJ);

opts.units = 'p';
planet3D('Earth', opts, rdE);
planet3D('Earth', opts, rfbE);
planet3D('Earth', opts, raE);

planet3D( 'Venus', opts, rdV);
planet3D( 'Venus', opts, rfbV);
planet3D( 'Venus', opts, raV);

legend('Jupiter orbit', 'Transfer arc 1', 'Earth orbit', 'Transfer arc 2', 'Venus orbit');
text(rdJ(1)+2e7, rdJ(2)+2e7, rdJ(3)+2e7, '$J_{dep}$', 'color','k','FontSize',14);
text(rfbJ(1)+2e7, rfbJ(2)+2e7, rfbJ(3)+2e7, '$J_{fb}$', 'color','k','FontSize',14);
text(raJ(1)+2e7, raJ(2)+2e7, raJ(3)+2e7, '$J_{arr}$', 'color','k','FontSize',14);

text(rdE(1)+2e7, rdE(2)+2e7, rdE(3)+2e7, '$E_{dep}$', 'color','k','FontSize',14);
text(rfbE(1)+2e7, rfbE(2)+2e7, rfbE(3)+2e7, '$E_{fb}$', 'color','k','FontSize',14);
text(raE(1)+2e7, raE(2)+2e7, raE(3)+2e7, '$E_{arr}$', 'color','k','FontSize',14);

text(rdV(1)+2e7, rdV(2)+2e7, rdV(3)+2e7, '$V_{dep}$', 'color','k','FontSize',14);
text(rfbV(1)+2e7, rfbV(2)+2e7, rfbV(3)+2e7, '$V_{fb}$', 'color','k','FontSize',14);
text(raV(1)+2e7, raV(2)+2e7, raV(3)+2e7, '$V_{arr}$', 'color','k','FontSize',14);

xlabel("$r_x[km]$", 'FontSize', 18);
ylabel("$r_y[km]$", 'FontSize', 18);
zlabel("$r_z[km]$", 'FontSize', 18);
title('\textbf{Trajectory}', 'FontSize', 18);
grid on
axis equal

%% PORK-CHOP PLOTS

run porkchopplots.m

%% EARTH'S FLYBY

interplanetary = figure();
interplanetary.WindowState = 'maximized';
hold on
% outgoing hyperbola
plot3(stateH1(1:50,1), stateH1(1:50,2), stateH1(1:50,3), 'red', 'linewidth', 2);   
% incoming hyperbola
plot3(stateH2(1:50,1), stateH2(1:50,2), stateH2(1:50,3), 'blue', 'linewidth', 2); 
opts.units = 'km';
grid on
axis equal
% showing the perigee passage
scatter3(stateH1(1,1), stateH1(1,2),stateH1(1,3), 'magenta' ,'filled', 'LineWidth', 10);
% add arrow in the direction of the planet's velocity 
arrow3d([0 vfbE(1)*10^3],[0 vfbE(2)*10^3],[0 vfbE(3)*10^3], 9/10, 0.2*10^3, 0.8*10^3, 'k');
planet3D('Earth', opts, [0 0 0]);
title('\textbf{Flyby in Earth-centred frame parallel to HECI}', 'FontSize', 18)
xlabel('X [km]', 'FontSize', 18)
ylabel('Y [km]', 'FontSize', 18)
zlabel('Z [km]', 'FontSize', 18)
legend('incoming hyperbola', 'outgoing hyperbola', 'pericenter' , 'planet velocity')

%% GA and FMINCON CONVERGENCE

% if the methods chosen is tho one that, after a manual mesh, looks for the
% global minimum around the minimum of the discretized mesh, the
% convergengence graph can be plotted
if options.method == 1

    minimization = figure();
    minimization.WindowState = 'maximized';
    plot(1:optim.noptim, results.ga_Dv, 'b', 'LineWidth', 2)
    hold on
    plot(1:optim.noptim, results.fmincon_Dv, 'r', 'LineWidth', 2)
    grid on
    xlabel('optimization iteration [-]', 'FontSize', 18)
    ylabel('$\Delta v_{tot}$ [km/s]', 'FontSize', 18)
    title('\textbf{Global minimum iterative search convergence}', 'FontSize', 18)
    legend('Genetic Algorithm', 'Gradient Method')

end

%% 3D PORK CHOP PLOT

pcp3d = figure();
scatter3(temp.pcp3d_x, temp.pcp3d_y, temp.pcp3d_z, 80, temp.pcp3d_Dv, 'LineWidth', 2)
hold on
scatter3(temp.pcp3dNA_x, temp.pcp3dNA_y, temp.pcp3dNA_z, 50, temp.pcp3dNA_Dv, 'LineWidth', 1.5, 'Marker', '+')
datetick( 'x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
datetick( 'y', 'yyyy mmm dd', 'keeplimits');
datetick( 'z', 'yyyy mmm dd', 'keeplimits');
xtickangle(45);
Color = colorbar;
colormap('turbo')
xlabel('Jupiter departure date [MJD2000]', 'FontSize', 18);
ylabel('Earth flyby date [MJD2000]', 'FontSize', 18);
zlabel('Venus arrival date [MJD2000]', 'FontSize', 18)
zlim([datenum(2025,08,01), datenum(2065, 08, 01)]);
caxis([12 100])
set(Color,'YtickLabel');
Color.Label.String = '$\Delta v [km/s]$';
Color.Label.Interpreter = 'latex';
grid on
title('\textbf{Dates 3D space Pork Chop Plot}', 'FontSize', 18);

pcp3d = figure();
scatter3(temp.pcp3d_x, temp.pcp3d_y - temp.pcp3d_x, temp.pcp3d_z - temp.pcp3d_y, 80, temp.pcp3d_Dv, 'LineWidth', 2)
hold on
scatter3(temp.pcp3dNA_x, temp.pcp3dNA_y - temp.pcp3dNA_x, temp.pcp3dNA_z - temp.pcp3dNA_y, 50, temp.pcp3dNA_Dv, 'LineWidth', 1.5, 'Marker', '+')
datetick( 'x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
xtickangle(45);
Color = colorbar;
colormap('turbo')
xlabel('Jupiter departure date [MJD2000]', 'FontSize', 18);
ylabel('Jpiter - Earth ToF [days]', 'FontSize', 18);
zlabel('Earth - Venus ToF [days]', 'FontSize', 18);
xlim([datenum(2025,08,01), datenum(2065, 08, 01)]);
zlim([0, arrival.tof_max + 100]);
ylim([0, departure.tof_max + 100]);
caxis([12 100])
set(Color,'YtickLabel');
Color.Label.String = '$\Delta v [km/s]$';
Color.Label.Interpreter = 'latex';
grid on
title('\textbf{Dates 3D space Pork Chop Plot}', 'FontSize', 18);


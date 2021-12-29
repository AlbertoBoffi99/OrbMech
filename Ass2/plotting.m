% PLOTTING FILE
% In this file all graphs are plotted

%-------------------------------------------------------------------------%

%% NOMINAL ORBIT

% plotting the nominal orbit
fig.nominal = figure();
earth_sphere
hold on 
plot3(nominal.car(:,1), nominal.car(:,2), nominal.car(:,3), 'r', 'LineWidth', 2 )
xlabel('$r_x [km]$', 'FontSize', 18)
ylabel('$r_y [km]$', 'FontSize', 18)
zlabel('$r_z [km]$', 'FontSize', 18)
title ('\textbf{Nominal orbit in an unperturbed condition}', 'FontSize', 18) 
axis equal  
grid on

%% GROUNDTRACKS

% first GT plot
fig.GT1 = figure();
image_file = 'earth_2D.jpg'; 
cdata = imread(image_file);
image([-180 180], [90 -90], cdata);
ax = gca;
ax.YDir = 'normal';

hold on 

% plot of unperturbed case
plot(wrapTo180(rad2deg(nominal.lon)), rad2deg(nominal.lat),'.', 'Color','g');
plot(wrapTo180(rad2deg(nominal.lon(1))), rad2deg(nominal.lat(1)),'*', 'LineWidth', 10, 'Color', 'g');
plot(wrapTo180(rad2deg(nominal.lon(end))), rad2deg(nominal.lat(end)),'^', 'LineWidth', 6, 'Color', 'g');

% plot of repeating unperturbed case
plot(wrapTo180(rad2deg(repun.lon)), rad2deg(repun.lat),'.', 'Color','r');
plot(wrapTo180(rad2deg(repun.lon(1))), rad2deg(repun.lat(1)),'*', 'LineWidth', 10, 'Color', 'r');
plot(wrapTo180(rad2deg(repun.lon(end))), rad2deg(repun.lat(end)),'^', 'LineWidth', 6, 'Color', 'r');

% plot of perturbed case (J2 + SRP)
plot(wrapTo180(rad2deg(per.lon)), rad2deg(per.lat),'.', 'Color','b');
plot(wrapTo180(rad2deg(per.lon(1))), rad2deg(per.lat(1)),'*', 'LineWidth', 10, 'Color', 'b');
plot(wrapTo180(rad2deg(per.lon(end))), rad2deg(per.lat(end)),'^', 'LineWidth', 6, 'Color', 'b');

legend('Unperturbed', 'START', 'END', ...
    'Repeating Unperturbed', 'START', 'END',...
    'Perturbed', 'START', 'END');

axis ([-180 180 -90 90]);
axis equal;
grid on;
xlabel('lon [deg]');
ylabel('lat [deg]');

%-------------------------------------------------------------------------%

% second GT plot
fig.GT2 = figure();

% to plot the earth as background
image_file = 'earth_2D.jpg'; % saved in the same folder
cdata = imread(image_file);
image([-180 180],[90 -90], cdata);
ax = gca;
ax.YDir = 'normal'; % to rotate the earth draw

hold on

% plot of repeating unperturbed case
plot(wrapTo180(rad2deg(repun.lon)), rad2deg(repun.lat),'.', 'Color','r');
plot(wrapTo180(rad2deg(repun.lon(1))), rad2deg(repun.lat(1)),'*', 'LineWidth', 10, 'Color', 'r');
plot(wrapTo180(rad2deg(repun.lon(end))), rad2deg(repun.lat(end)),'^', 'LineWidth', 6, 'Color', 'r');

% plot of repeating perturbed case
plot(wrapTo180(rad2deg(repper.lon)), rad2deg(repper.lat),'.', 'Color','y');
plot(wrapTo180(rad2deg(repper.lon(1))), rad2deg(repper.lat(1)),'*', 'LineWidth', 10, 'Color', 'y');
plot(wrapTo180(rad2deg(repper.lon(end))), rad2deg(repper.lat(end)),'^', 'LineWidth', 6, 'Color', 'y');

legend('Repeating Unperturbed', 'START', 'END', ...
    'Repeating Perturbed', 'START', 'END');

axis ([-180 180 -90 90]);
axis equal;
grid on;
xlabel('lon [deg]');
ylabel('lat [deg]');

%% PERTURBED ORBITS

% https://it.mathworks.com/matlabcentral/answers/101346-how-do-i-use-multiple-colormaps-in-a-single-figure-in-r2014a-and-earlier#Example_1

% propagation of the orbit in Cartesian elements
fig.per = figure();
patch(per.car(:,1), per.car(:,2), per.car(:,3), per.time_car, 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1);
c = colorbar('TickLabelInterpreter', 'latex');
c.Label.String = '[s]';
c.Label.FontSize = 14;
c.Label.Interpreter = 'latex';
hold off
opts.units = 'km';
planet3D('Earth', opts, [0 0 0]);
xlabel('$r_x [km]$')
ylabel('$r_y [km]$')
zlabel('$r_z [km]$')
axis equal
grid on
title("Propagation of the perturbed orbit in Cartesian elements");

% propagation of the orbit in Cartesian elements converted from gauss
fig.per_gauss = figure();
patch(per.car_gauss(:,1), per.car_gauss(:,2), per.car_gauss(:,3), per.time_kep, 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1);
c = colorbar('TickLabelInterpreter', 'latex');
c.Label.String = '[s]';
c.Label.FontSize = 14;
c.Label.Interpreter = 'latex';
hold off
opts.units = 'km';
planet3D('Earth', opts, [0 0 0]);
xlabel('$r_x [km]$')
ylabel('$r_y [km]$')
zlabel('$r_z [km]$')
axis equal
grid on
title("Propagation of the perturbed orbit in Cartesian elements from Keplerian elements");

%-------------------------------------------------------------------------%

% relative errors plot

for i = [1 2]

    switch i
        case 1
            name = 'a';
            pos = 1;
            denom = nominal.kep(1);
        
        case 2
            name = 'e';
            pos = 2;
            denom = 1;
            
    end
    
    figure();
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, ...
        abs(per.kep_unwr(:,pos) - per.kep_gauss(:,pos))/denom);
    grid on
    % title ("relative error of 'name' between gauss and cartesian propagation")
    
    figure();
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, per.kep_gauss(:,pos));
    hold on
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, per.kep_unwr(:,pos));
    grid on
    legend ('Gauss', 'Cartesian')
    titlestr = strcat ("Propagation of ", name);
    title(titlestr);

end

for i = 1 : 4

    switch i
         case 1
            name = 'i';
            pos =3;
            denom = 2*pi;
    
        case 2
            name = 'Om';
            pos = 4;
            denom = 2*pi;
    
        case 3
            name = 'om';
            pos = 5;
            denom = 2*pi;
    
        case 4
            name = 'f';
            pos = 6;
            denom = per.kep_unwr(:,6);
    
    end    
    
    figure();
    plot (temp.tspan_3year/(24*3600)-real.mjd2000_start, rad2deg(abs(per.kep_unwr(:,pos) - per.kep_gauss(:,pos)))./abs(denom));
    grid on
    titlestr = strcat("Relative error of ", name, " between Gauss and cartesian propagation");
    title (titlestr)
    
    figure();
    plot (temp.tspan_3year/(24*3600)-real.mjd2000_start, rad2deg(per.kep_gauss(:,pos)));
    hold on
    plot (temp.tspan_3year/(24*3600)-real.mjd2000_start, rad2deg(per.kep_unwr(:,pos)));
    grid on
    legend ('Gauss', 'Cartesian')

end

%% GAUSS PERTURBED VS REAL TRAJECTORY

% condizioni iniziali nuova integrazione su tspan_real come quelle del
% satellite reale (f0 arbitraria)

% plotta confronto reale vs gauss perturbed

for i = 1:2
    gaussvsreal = figure();
    plot (temp.tspan_real/(24*3600)-real.mjd2000_start, per.kep_gauss_real(:,i), ...
          temp.tspan_real/(24*3600)-real.mjd2000_start, real.kep(:,i));
    grid on
    legend ('Gauss', 'Real')
end

for i = 3:5
    gaussvsreal = figure();
    plot (temp.tspan_real/(24*3600)-real.mjd2000_start, rad2deg(per.kep_gauss_real(:,i)), ...
        temp.tspan_real/(24*3600)-real.mjd2000_start, real.kep(:,i));
    grid on
    legend ('Gauss', 'Real')
end


%% FILTERING OF PERTURBED KEPLERIAN ELEMETS

for i = 1 : 2

    switch i
        case 1
            pos = 1;
            % number of elements in the movemean
            filter_par = nominal.T/(3600*24*365*3)*length(temp.tspan_3year);
        
        case 2
            pos = 2;
            % number of elements in the movemean
            filter_par = nominal.T/(3600*24*365*3)*length(temp.tspan_3year);
    end
    
    filt.kep_gauss = movmean(per.kep_gauss(:,pos), filter_par);
    
    figure();
    plot (temp.tspan_3year/(24*3600)-real.mjd2000_start, per.kep_gauss(:,pos));
    hold on
    plot (temp.tspan_3year/(24*3600)-real.mjd2000_start, filt.kep_gauss);
    grid on
    legend ('Original', 'Filtered')

end

for i = 3 : 6

    switch i
        case 3
            name = 'i';
            pos = 3;
            % number of elements in the movemean
            filter_par = nominal.T/(3600*24*365*3)*length(temp.tspan_3year);
        
        case 4
            name = 'Om';
            pos = 4;
            % number of elements in the movemean
            filter_par = nominal.T/(3600*24*365*3)*length(temp.tspan_3year);
    
        case 5
            name = 'om';
            pos = 5;
            % number of elements in the movemean
            filter_par = nominal.T/(3600*24*365*3)*length(temp.tspan_3year);
    
        case 6
            name = 'f';
            pos = 6;
            % number of elements in the movemean
            filter_par = nominal.T/(3600*24*365*3)*length(temp.tspan_3year);
    end
    
    filt.kep_gauss = movmean(per.kep_gauss(:,pos), filter_par);

    figure();
    plot (temp.tspan_3year/(24*3600)-real.mjd2000_start, rad2deg(per.kep_gauss(:,pos)));
    hold on
    plot (temp.tspan_3year/(24*3600)-real.mjd2000_start, rad2deg(filt.kep_gauss));
    grid on
    legend ('Original', 'Filtered')

end
% PLOTTING FILE
% In this file all graphs are plotted

%-------------------------------------------------------------------------%

%% COLORS

% dark red
nominal.color = [172 5 5]/256;
% light orange
repun.color = [255 142 27]/256;
% light red
repper.color = [255 55 89]/256;
% dark light blu ;)
per.color = [0 178 172]/256;
% light purple
real.color = [204 134 255]/256;
% brown
filt.color = [81 53 26]/256;

%% NOMINAL ORBIT

% plotting the nominal orbit
figure();
set(gcf, 'WindowState', 'Maximized')
opts.units = 'km';
planet3D('Earth', opts, [0 0 0]);
hold on 
plot3(nominal.car(:,1), nominal.car(:,2), nominal.car(:,3), 'Color', nominal.color, 'LineWidth', 2 )
xlabel('$r_x [km]$', 'FontSize', 14)
ylabel('$r_y [km]$', 'FontSize', 14)
zlabel('$r_z [km]$', 'FontSize', 14)
title ('\textbf{Nominal orbit in an unperturbed condition}', 'FontSize', 18)
subtitle_str = strcat("$a$ = ", num2str(nominal.a0), " km;",...
                  "\hspace{1cm}$e$ = ",  num2str(nominal.e0), ...
                  "\hspace{1cm}$i$ = ",  num2str(nominal.i0), " deg;", ...
                  "\hspace{1cm}$\Omega$ = ",  num2str(nominal.Ome0), " deg;", ...
                  "\hspace{1cm}$\omega$ = ",  num2str(nominal.ome0), " deg;", ...
                  "\hspace{1cm}$\theta$ = ",  num2str(nominal.f0), " deg;");
subtitle(subtitle_str, 'FontSize', 14)
axis equal  
grid on

% save figure
if options.save
    savefig('.\Results\nominal');
end

%% GROUND TRACKS

figure();
set(gcf, 'WindowState', 'Maximized')
% plot of unperturbed case
hold off
geoplot( rad2deg(nominal.lat), wrapTo180(rad2deg(nominal.lon)), 'o', 'MarkerSize', 2, 'Color', nominal.color);
geolimits([-50 50], [-180 180])
geobasemap darkwater
hold on
set(gca,'FontSize',14)
set(gca,'FontName', 'Computer Modern')
geoplot(rad2deg(nominal.lat(1)), wrapTo180(rad2deg(nominal.lon(1))), '*', 'MarkerSize', 20, 'LineWidth', 2, 'Color', nominal.color);
geoplot(rad2deg(nominal.lat(end)), wrapTo180(rad2deg(nominal.lon(end))), 'o', 'MarkerSize', 12, 'LineWidth', 2, 'Color', nominal.color);

% plot of repeating unperturbed case
geoplot( rad2deg(repun.lat), wrapTo180(rad2deg(repun.lon)), 'o', 'MarkerSize', 2, 'Color', repun.color);
geoplot(rad2deg(repun.lat(1)), wrapTo180(rad2deg(repun.lon(1))), '*', 'MarkerSize', 20, 'LineWidth', 2, 'Color', repun.color);
geoplot(rad2deg(repun.lat(end)), wrapTo180(rad2deg(repun.lon(end))), 'o', 'MarkerSize', 12, 'LineWidth', 2, 'Color', repun.color);

% plot of perturbed case (J2 + SRP)
geoplot(rad2deg(per.lat), wrapTo180(rad2deg(per.lon)), 'o', 'MarkerSize', 2, 'Color', per.color);
geoplot(rad2deg(per.lat(1)), wrapTo180(rad2deg(per.lon(1))), '*', 'MarkerSize', 20, 'LineWidth', 2, 'Color', per.color);
geoplot(rad2deg(per.lat(end)), wrapTo180(rad2deg(per.lon(end))), 'o', 'MarkerSize', 12, 'LineWidth', 2, 'Color', per.color);

lgd = legend('Unperturbed', 'Start', 'End', ...
    'Repeating Unperturbed', 'Start', 'End',...
    'Perturbed', 'Start', 'End');
lgd.Location = 'southeastoutside';

title ('\textbf{Unperturbed, unperturbed and repeating, and perturbed ground tracks}', 'FontSize', 18)
subtitle_str = strcat('\textit{Time period: }', temp.periodstr);
subtitle(subtitle_str, 'FontSize', 14)

% save figure
if options.save
    savefig('.\Results\GT1');
end

%-------------------------------------------------------------------------%

figure();
set(gcf, 'WindowState', 'Maximized')
% plot of repeating unperturbed case
hold off
geoplot(rad2deg(repper.lat), wrapTo180(rad2deg(repper.lon)), 'o', 'MarkerSize', 2, 'Color', repper.color);
geolimits([-50 50], [-180 180])
geobasemap darkwater
hold on
set(gca,'FontSize',14)
set(gca,'FontName', 'Comclcuter Modern')
geoplot(rad2deg(repper.lat(1)), wrapTo180(rad2deg(repper.lon(1))), '*', 'MarkerSize', 20, 'LineWidth', 2, 'Color', repper.color);
geoplot(rad2deg(repper.lat(end)), wrapTo180(rad2deg(repper.lon(end))), 'o', 'MarkerSize', 12, 'LineWidth', 2, 'Color', repper.color);
geoplot(rad2deg(repun.lat(1)), wrapTo180(rad2deg(repun.lon(1))), '*', 'MarkerSize', 20, 'LineWidth', 2, 'Color', repun.color);
geoplot(rad2deg(repun.lat(end)), wrapTo180(rad2deg(repun.lon(end))), 'o', 'MarkerSize', 12, 'LineWidth', 2, 'Color', repun.color);
% plot of repeating perturbed case
geoplot( rad2deg(repun.lat), wrapTo180(rad2deg(repun.lon)), 'o', 'MarkerSize', 2, 'Color', repun.color);
geoplot(rad2deg(repun.lat(1)), wrapTo180(rad2deg(repun.lon(1))), '*', 'MarkerSize', 20, 'LineWidth', 2, 'Color', repun.color);
geoplot(rad2deg(repun.lat(end)), wrapTo180(rad2deg(repun.lon(end))), 'o', 'MarkerSize', 12, 'LineWidth', 2, 'Color', repun.color);
lgd = legend('Repeating Perturbed', 'Start', 'End', ...
    'Repeating Unperturbed', 'Start', 'End');
lgd.Location = 'southeastoutside';
lgd.FontSize = 16;

title ('\textbf{Repeating perturbed and uperturbed ground tracks}', 'FontSize', 18)
subtitle_str = strcat('\textit{Time period: }', temp.periodstr);
subtitle(subtitle_str, 'FontSize', 14)

% save figure
if options.save
    savefig('.\Results\GT2');
end

%% PERTURBED ORBITS

% https://it.mathworks.com/matlabcentral/answers/101346-how-do-i-use-multiple-colormaps-in-a-single-figure-in-r2014a-and-earlier#Example_1

% propagation of the orbit in Cartesian elements
figure();
set(gcf, 'WindowState', 'Maximized')
p = patch(per.car3D(:,1), per.car3D(:,2), per.car3D(:,3), per.time_car3D(:)/(24*3600), 'FaceColor', 'none', 'EdgeColor', 'interp');
c = colorbar('TickLabelInterpreter', 'latex');
colormap('turbo')
c.Label.String = '[MJD2000]';
c.Label.FontSize = 14;
c.Label.Interpreter = 'latex';
set(c.XLabel,{'Rotation','Position'},{0,[0.5 -0.1]})
set(gca, 'FontSize', 14)
hold off
opts.units = 'km';
planet3D('Earth', opts, [0 0 0]);
xlabel('$r_x [km]$', 'FontSize', 14)
ylabel('$r_y [km]$', 'FontSize', 14)
zlabel('$r_z [km]$', 'FontSize', 14)
axis equal
grid on
title("\textbf{Propagation of the J2 and SRP perturbed orbit}", 'FontSize', 18);
subtitle("\textit{Integration performed in Cartesian elements}", 'FontSize', 14);

% save figure
if options.save
    savefig('.\Results\propcar');
end

% Gauss propagation of the orbit elements
figure();
set(gcf, 'WindowState', 'Maximized')
patch(per.car_gauss3D(:,1), per.car_gauss3D(:,2), per.car_gauss3D(:,3), per.time_kep3D/(24*3600), 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1);
c = colorbar('TickLabelInterpreter', 'latex');
colormap('turbo')
c.Label.String = '[MJD2000]';
c.Label.FontSize = 14;
c.Label.Interpreter = 'latex';
set(c.XLabel,{'Rotation','Position'},{0,[0.5 -0.1]})
set(gca, 'FontSize', 14)
hold off
opts.units = 'km';
planet3D('Earth', opts, [0 0 0]);
xlabel('$r_x [km]$', 'FontSize', 14)
ylabel('$r_y [km]$', 'FontSize', 14)
zlabel('$r_z [km]$', 'FontSize', 14)
axis equal
grid on
title("\textbf{Propagation of the J2 and SRP perturbed orbit}", 'FontSize', 18);
subtitle("\textit{Integration performed through Gauss propagation of Keplerian elements}", 'FontSize', 14);

% save figure
if options.save
    savefig('.\Results\propGauss');
end

%-------------------------------------------------------------------------%

%% RELATIVE ERRORS

for i = [1 2]

    switch i
        case 1
            titlestr = "\textbf{Propagation of semi-major axis: $a$}";
            subtitlestr = "\textit{semi-major axis: } $a$";
            ylabelstr = "a [km]";
            denom = nominal.kep(1);
        
        case 2
            titlestr = "\textbf{Propagation of eccentricity: $e$}";
            subtitlestr = "\textit{eccentricity: } $e$";
            ylabelstr = "e [-]";
            denom = 1;
            
    end
    
    figure();
    set(gcf, 'WindowState', 'Maximized')
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, ...
        abs(per.kep_unwr(:,i) - per.kep_gauss(:,i))*100/denom, 'k', 'LineWidth', 2 );
    grid on
    set(gca, 'FontSize', 14);
    xlabel('time [MJD2000]', 'FontSize', 14)
    ylabel("error \%", 'FontSize', 14)
    title("\textbf{Relative error between Gauss and Cartesian propagation}", 'FontSize', 18);
    subtitle(subtitlestr, "FontSize", 14)

    % save figure
    if options.save
        filename = strcat('.\Results\err_rel', num2str(i));
        savefig(filename);
    end
    
    figure();
    set(gcf, 'WindowState', 'Maximized')
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, per.kep_gauss(:,i), 'Color', repper.color, 'LineWidth', 2 );
    hold on
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, per.kep_unwr(:,i), 'Color', repun.color, 'LineWidth', 2 );
    grid on
    set(gca, 'FontSize', 14);
    xlabel('time [MJD2000]', 'FontSize', 14)
    ylabel(ylabelstr, 'FontSize', 14)
    lgd = legend ('Gauss', 'Cartesian');
    title(titlestr, 'FontSize', 18);
    lgd.Location = 'southeastoutside';
    lgd.FontSize = 16;

    % save figure
    if options.save
        filename = strcat('.\Results\el_prop', num2str(i));
        savefig(filename);
    end

end

for i = 3 : 6

    switch i
         case 3
            titlestr = "\textbf{Propagation of inclination: $i$}";
            subtitlestr = "\textit{inclination: } $i$";
            ylabelstr = "i [deg]";
            denom = 2*pi;
    
        case 4
            titlestr = "\textbf{Propagation of RAAN: $\Omega$}";
            subtitlestr = "\textit{RAAN: } $\Omega$";
            ylabelstr = "$\Omega$ [deg]";
            denom = 2*pi;
    
        case 5
            titlestr = "\textbf{Propagation of pericenter anomaly: $\omega$}";
            subtitlestr = "\textit{pericenter anomaly: } $\omega$";
            ylabelstr = "$\omega$ [deg]";
            denom = 2*pi;
    
        case 6
            titlestr = "\textbf{Propagation of true anomaly: $\theta$}";
            subtitlestr = "\textit{true anomaly: } $\theta$";
            ylabelstr = "$\theta$ [deg]";
            denom = per.kep_unwr(:,6);
    
    end    
    
    figure();
    set(gcf, 'WindowState', 'Maximized')
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, ...
        abs(per.kep_unwr(:,i) - per.kep_gauss(:,i))*100./denom, 'k', 'LineWidth', 2 );
    grid on
    set(gca, 'FontSize', 14);
    xlabel('time [MJD2000]', 'FontSize', 14)
    ylabel("error \%", 'FontSize', 14)
    title("\textbf{Relative error between Gauss and Cartesian propagation}", 'FontSize', 18);
    subtitle(subtitlestr, "FontSize", 14)

%     % save figure
%     if options.save
%         filename = strcat('.\Results\err_rel', num2str(i));
%         savefig(filename);
%     end
    
    figure();
    set(gcf, 'WindowState', 'Maximized')
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, rad2deg(per.kep_gauss(:,i)), 'Color', repper.color, 'LineWidth', 2 );
    hold on
    plot(temp.tspan_3year/(24*3600)-real.mjd2000_start, rad2deg(per.kep_unwr(:,i)), 'Color', repun.color, 'LineWidth', 2 );
    grid on
    set(gca, 'FontSize', 14);
    xlabel('time [MJD2000]', 'FontSize', 14)
    ylabel(ylabelstr, 'FontSize', 14)
    lgd = legend ('Gauss', 'Cartesian');
    title(titlestr, 'FontSize', 18);
    lgd.Location = 'southeastoutside';
    lgd.FontSize = 16;

%     % save figure
%     if options.save
%         filename = strcat('.\Results\el_prop', num2str(i));
%         savefig(filename);
%     end

end

%% GAUSS PERTURBED VS REAL TRAJECTORY

% plotta confronto reale vs gauss perturbed

figure();
set(gcf, 'WindowState', 'Maximized')
plot (per.time_kep_real/(24*3600)-real.mjd2000_start, per.kep_gauss_real_filt(:,1), 'Color', repper.color, 'LineWidth', 2 );
hold on
plot (temp.tspan_real/(24*3600)-real.mjd2000_start, real.kep(:,1), 'Color', real.color, 'LineWidth', 2 );
grid on
lgd = legend ('Gauss', 'Real');
set(gca, 'FontSize', 14);
xlabel('time [MJD2000]', 'FontSize', 14)
ylabel("a [km]", 'FontSize', 14)
title("\textbf{Comparison between Gauss-propagated and real perturbated orbits}", 'FontSize', 18);
subtitle("\textit{semi-major axis: } $a$", "FontSize", 14)
lgd.Location = 'southeastoutside';
lgd.FontSize = 16;

% save figure
if options.save
    savefig('.\Results\gauss_vs_real1');
end

figure();
set(gcf, 'WindowState', 'Maximized')
plot (per.time_kep_real/(24*3600)-real.mjd2000_start, per.kep_gauss_real(:,2), 'Color', repper.color, 'LineWidth', 2 );
hold on
plot (temp.tspan_real/(24*3600)-real.mjd2000_start, real.kep(:,2), 'Color', real.color, 'LineWidth', 2 );
grid on
lgd = legend ('Gauss', 'Real');
set(gca, 'FontSize', 14);
xlabel('time [MJD2000]', 'FontSize', 14)
ylabel("e [-]", 'FontSize', 14)
title("\textbf{Comparison between Gauss-propagated and real perturbated keplerian elements}", 'FontSize', 18);
subtitle("\textit{eccentricity: } $e$", "FontSize", 14)
lgd.Location = 'southeastoutside';
lgd.FontSize = 16;

% save figure
if options.save
    savefig('.\Results\gauss_vs_real2');
end

for i = 3:5

    switch i
         case 1
            subtitlestr = "\textit{inclination: } $i$";
            ylabelstr = "i [deg]";
    
        case 2
            subtitlestr = "\textit{RAAN: } $\Omega$";
            ylabelstr = "$\Omega$ [deg]";
    
        case 3
            subtitlestr = "\textit{pericenter anomaly: } $\omega$";
            ylabelstr = "$\omega$ [deg]";
    end

    figure();
    set(gcf, 'WindowState', 'Maximized')
    plot (per.time_kep_real/(24*3600)-real.mjd2000_start, rad2deg(per.kep_gauss_real(:,i)), 'Color', repper.color, 'LineWidth', 2 );
    hold on
    plot (temp.tspan_real/(24*3600)-real.mjd2000_start, real.kep(:,i), 'Color', real.color, 'LineWidth', 2 );
    grid on
    lgd = legend ('Gauss', 'Real');
    set(gca, 'FontSize', 14);
    xlabel('time [MJD2000]', 'FontSize', 14)
    ylabel(ylabelstr, 'FontSize', 14)
    title("\textbf{Comparison between Gauss-propagated and real perturbated keplerian elements}", 'FontSize', 18);
    subtitle(subtitlestr, "FontSize", 14)
    lgd.Location = 'southeastoutside';
    lgd.FontSize = 16;

    % save figure
    if options.save
        filename = strcat('.\Results\gauss_vs_real', num2str(i));
        savefig(filename);
    end
end


%% FILTERING OF PERTURBED KEPLERIAN ELEMENTS

HFfilter = nominal.T/(5*24*3600)*length(per.time_kep_short);

for i = 1 : 2

    switch i
        case 1
            subtitlestr = "\textit{semi-major axis: } $a$";
            ylabelstr = "a [km]";
        
        case 2
            subtitlestr = "\textit{eccentricity: } $e$";
            ylabelstr = "e [-]";
            denom = 1;
            
    end
    
    filt.kep_gauss = movmean(per.kep_gauss_short(:,i), HFfilter);
    
    figure();
    set(gcf, 'WindowState', 'Maximized')
    plot (temp.tspan_5day(0.5/5*40000:end-0.5/5*40000)/(24*3600)-real.mjd2000_start, ...
        per.kep_gauss_short(0.5/5*40000:end-0.5/5*40000,i), 'Color', repper.color, 'LineWidth', 2 );
    hold on
    plot (temp.tspan_5day(0.5/5*40000:end-0.5/5*40000)/(24*3600)-real.mjd2000_start, ...
        filt.kep_gauss(0.5/5*40000:end-0.5/5*40000), 'Color', filt.color, 'LineWidth', 2 );
    grid on
    lgd = legend ('Original Gauss', 'Filtered');
    lgd.Location = 'southeastoutside';
    lgd.FontSize = 16;
    set(gca, 'FontSize', 14);
    xlabel('time [MJD2000]', 'FontSize', 14)
    ylabel(ylabelstr, 'FontSize', 14)
    title("\textbf{Comparison between Gauss-propagated and filtered keplerian elements}", 'FontSize', 18);
    subtitle(subtitlestr, "FontSize", 14)
    xlim([0.5 4.5])

    % save figure
    if options.save
        filename = strcat('.\Results\filter', num2str(i));
        savefig(filename);
    end

end

for i = 3 : 6

    switch i
         case 3
            subtitlestr = "\textit{inclination: } $i$";
            ylabelstr = "i [deg]";
    
        case 4
            subtitlestr = "\textit{RAAN: } $\Omega$";
            ylabelstr = "$\Omega$ [deg]";
    
        case 5
            subtitlestr = "\textit{pericenter anomaly: } $\omega$";
            ylabelstr = "$\omega$ [deg]";

        case 6
            subtitlestr = "\textit{true anomaly: } $\theta$";
            ylabelstr = "$\theta$ [deg]";

    end
    
    filt.kep_gauss = movmean(per.kep_gauss_short(:,i), HFfilter);

    figure();
    set(gcf, 'WindowState', 'Maximized')
    plot (temp.tspan_5day(0.5/5*40000:end-0.5/5*40000)/(24*3600)-real.mjd2000_start, ...
        rad2deg(per.kep_gauss_short(0.5/5*40000:end-0.5/5*40000,i)), 'Color', repper.color, 'LineWidth', 2 );
    hold on
    plot (temp.tspan_5day(0.5/5*40000:end-0.5/5*40000)/(24*3600)-real.mjd2000_start, ...
        rad2deg(filt.kep_gauss(0.5/5*40000:end-0.5/5*40000)), 'Color', filt.color, 'LineWidth', 2 );
    grid on
    lgd = legend ('Original Gauss', 'Filtered');
    lgd.Location = 'southeastoutside';
    lgd.FontSize = 16;
    set(gca, 'FontSize', 14);
    xlabel('time [MJD2000]', 'FontSize', 14)
    ylabel(ylabelstr, 'FontSize', 14)
    title("\textbf{Comparison between Gauss-propagated and filtered keplerian elements}", 'FontSize', 18);
    subtitle(subtitlestr, "FontSize", 14)
    xlim([0.5 4.5])

    % save figure
    if options.save
        filename = strcat('.\Results\filter', num2str(i));
        savefig(filename);
    end
end
% ASSIGNMENT 1
%   minimize the delta velocity required to depart from Jupiter and arrive
%   on Venus through a fly-by around the Earth

%-------------------------------------------------------------------------%

% initialize MATLAB space
clear all
close all
clc

% clock started
tic
fprintf('Execution started ...\n');

%-------------------------------------------------------------------------%

%% INPUT

% calling global functions
global out

% optimization method variable
options.method = 2;

% run configuration file, containing all numerical data
run config.m

%%  OPTIMIZATION

% declare flyby modelling function as handle
JEV = @(dates) JEV_lamb_fb_lamb(dates, astro, options);

% discretize departure window
departure.window = linspace(departure.date_min, arrival.date_max, optim.nmesh_dep);

% open waitbar
temp.wb = waitbar(0, 'Optimization Running ...');

% for every departure date inside the departure window
for j = 1:1:size(departure.window,2)-2;
    
    % updating waitbar
    waitbar(j/size(departure.window,2), temp.wb);

   % B constraint matrix for GA and fmincon
    optim.Bcon = [ - departure.window(j); ...
           - departure.window(j) - departure.tof_min; ...
           - departure.window(j) - departure.tof_min - arrival.tof_min; ...
                    departure.window(j+1); ...
             departure.window(j+1) + departure.tof_max; ...
             departure.window(j+1) + departure.tof_max + arrival.tof_max; ...
            - departure.date_min; ...
            arrival.date_max;
            optim.tof_minimum;
            optim.tof_minimum];

    % optimization of Dv using fmincon
    [temp.dates, temp.Dv, temp.exitflag] = ...
    ga(JEV, optim.ga_nvars, optim.Acon, optim.Bcon, ...
    [], [], [], [], [], options.ga_options);

    results.ga_Dv(j) = temp.Dv;
    results.ga_dates(j, 1:3) = temp.dates;
    results.ga_exitflag(j) = temp.exitflag;
    
    % if the delta velocity found is the minimum until now, save it
    % as the best result
    if temp.Dv < results.Dv && ~out.nanflag
        results.dates = temp.dates;
        results.Dv = temp.Dv;
        results.ViT1 = out.ViT1;
        results.ViT2 = out.ViT2;
        results.DvT1 = out.DvT1;
        results.DvT2 = out.DvT2;
        results.DvGA = out.DvGA;
        results.Dvfb = out.Dvfb;
        results.rp = out.rp;
        results.rp_norm = out.rp_norm;
        results.hp = results.rp_norm - astro.RE;
        results.vp_p = out.vp_p;
        results.vp_m = out.vp_m;
    end
end

% close waitbar
delete(temp.wb)

fprintf('\n     Optimization with discretization mesh conlcuded\n')

%% TIME SPENT DURING FLYBY

run fbtime.m

%-------------------------------------------------------------------------%

% clock stopped
fprintf('\nExecution concluded\n')
toc

%% PLOTTING

% if the user asks for plotting 
if options.plot

    fprintf('\nPlotting started ...\n')

    % plotting function
    run plotting.m
    
    fprintf('Plotting concluded\n')

end

%% DATA RECORD

% convert dates into normal date format from mjd2000
temp.dates = results.dates;
results.dates(1,1:6) = mjd20002date(temp.dates(1));
results.dates(2,:) = mjd20002date(temp.dates(2));
results.dates(3,:) = mjd20002date(temp.dates(3));

fprintf('\nDATA RECORD:\n')

fprintf('     departure date: %g/%g/%g [mm/dd/yyyy] \n', results.dates(1,2), results.dates(1,3), results.dates(1,1));
fprintf('     fly-by date: %g/%g/%g [mm/dd/yyyy] \n', results.dates(2,2), results.dates(2,3), results.dates(2,1));
fprintf('     arrival date: %g/%g/%g [mm/dd/yyyy] \n\n', results.dates(3,2), results.dates(3,3), results.dates(3,1));

fprintf('     total delta velocity needed: %g [km/s] \n', results.Dv);
fprintf('     delta velocity to begin interplanetary trajectory: %g [km/s] \n', results.DvT1);
fprintf('     natural fly-by delta velocity: %g [km/s] \n', results.Dvfb);
fprintf('     powered fly-by delta velocity gven at the perigee: %g [km/s] \n', results.DvGA);
fprintf('     delta velocity to end interplanetary trajectory: %g [km/s] \n\n', results.DvT2);

fprintf('     height of the perigee of the fly-by trajectory: %g [km] \n', results.hp);
fprintf("     time spent inside Earth's SOI: %g [s] \n\n", results.Dtfb);

%% SAVE RESULTS

% if the user asks for saving results 
if options.save

    save('.\Results\meshwga_results.mat', 'results');
    savefig(flyby, '.\Figures\meshwga\flyby.fig');
    savefig(pcpL1, '.\Figures\meshwga\pcpL1.fig');
    savefig(pcpL2, '.\Figures\meshwga\pcpL2.fig');
    savefig(interplanetary, '.\Figures\meshwga\interplanetary.fig');
    fprintf('Results and figure saved\n');

end

% clear workspace
clearvars -except optim astro options results departure arrival fly



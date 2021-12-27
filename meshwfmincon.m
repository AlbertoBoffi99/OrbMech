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
for j = 1:1:size(departure.window,2);
    
    % define the flyby window according to ToF retreived from 1st arc's pork
    % chop plot
    fly.window = linspace(departure.window(j) + departure.tof_min, ...
        departure.window(j) + departure.tof_max, optim.nmesh_fly);
    
    % updating waitbar
    waitbar(j/size(departure.window,2), temp.wb);

    % for every flyby date inside the flyby window
    for k = 1:1:size(fly.window,2);
    
        % define the arrival window according to ToF retreived from 2nd arc's pork
        % chop plot
        arrival.window = linspace(fly.window(k) + arrival.tof_min, ...
            min(fly.window(k) + arrival.tof_max, arrival.date_max), optim.nmesh_arr);

        % for every arrival date insisde the arrival window
        for z = 1:1:size(arrival.window,2);
            
            % departure, flyby, arrival dates vector
            temp.dates = [departure.window(j), fly.window(k), arrival.window(z)];

           % B constraint matrix for GA and fmincon
            optim.Bcon = [ - results.dates(1) + departure.SW; ...
                   - results.dates(2) + fly.SW; ...
                   - results.dates(3) + arrival.SW; ...
                     results.dates(1) + departure.SW; ...
                     results.dates(2) + fly.SW; ...
                     results.dates(3) + arrival.SW; ...
                    - departure.date_min; ...
                    arrival.date_max;
                    0;
                    0];
    
            % optimization of Dv using fmincon
            [temp.dates, temp.Dv, temp.exitflag] = ...
            fmincon(JEV, temp.dates, optim.Acon, optim.Bcon, ...
            [], [], [], [], [], options.fmincon_options);

            temp.Dv
            temp.exitflag
            
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
    end
end

% close waitbar
delete(temp.wb)

fprintf('\n     Optimization with discretization mesh conlcuded\n')

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

    save('.\Results\meshwfmincon_results.mat', 'results');
    fprintf('Results saved\n');

end

% clear workspace
clearvars -except optim astro options results departure arrival fly



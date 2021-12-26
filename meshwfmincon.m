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

% run configuration file, containing all numerical data
run config.m

%%  MESH OPTIMIZATION

% declare flyby modelling function as handle
JEV = @(dates) JEV_lamb_fb_lamb(dates, astro, options);

% non linear function for GA and fmincon
optim.nonlincon = @(x) rp_nonlinFcn(x, astro.RE, astro.h_atm);

% discretize departure window
departure.window = linspace(departure.date_min, departure.date_max, optim.nmesh_dep);

% open waitbar
temp.wb = waitbar(0, 'Pattern Optimization Running ...');

% for every departure date inside the departure window
for j = 1:1:size(departure.window,2);
    
    % define the flyby window according to ToF retreived from 1st arc's pork
    % chop plot
    fly.window = linspace(departure.window(j) + departure.tof_min, ...
        departure.window(j) + departure.tof_max, optim.nmesh_fly);
    
    % update waitbar
    waitbar(j/size(departure.window,2), temp.wb);
    
    % for every flyby date inside the flyby window
    for k = 1:1:size(fly.window,2);
        
        % define the arrival window according to ToF retreived from 2nd arc's pork
        % chop plot
        arrival.window = linspace(fly.window(k) + arrival.tof_min, ...
            fly.window(k) + arrival.tof_max, optim.nmesh_arr);
        
        % for every arrival date insisde the arrival window
        for z = 1:1:size(arrival.window,2);
            
            % departure, flyby, arrival dates vector
            temp.dates = [departure.window(j), fly.window(k), arrival.window(z)];

            % B constraint matrix for GA and fmincon
            optim.Bcon = [ - temp.dates(1) + departure.SW; ...
            - temp.dates(2) + fly.SW; ...
            - temp.dates(3) + arrival.SW; ...
            temp.dates(1) + departure.SW; ...
            temp.dates(2) + fly.SW; ...
            temp.dates(3) + arrival.SW ];

            % update fmincon guess according to ga result
            fmincon_guess = temp.dates;
    
            % optimization of Dv using fmincon
            [temp.dates, temp.Dv] = ...
            fmincon(JEV, fmincon_guess, optim.Acon, optim.Bcon, ...
            [], [], [], [], optim.nonlincon, options.fmincon_options);
            
            % if the delta velocity found is the minimum until now, save it
            % as the best result
            if temp.Dv < results.Dv_min && ~out.nanflag
                results.dates = temp.dates;
                results.Dv_min = temp.Dv;
                results.ViT1 = out.ViT1;
                results.ViT2 = out.ViT2;
                results.DvT1 = out.DvT1;
                results.DvT2 = out.DvT2;
                results.Dvfb = out.Dvfb;
                results.rp_norm = out.rp_norm;
                results.rp = out.rp;
                results.vp_p = out.vp_p;
                results.vp_m = out.vp_m;
            end
        end
    end
end

% close waitbar
delete(temp.wb)
% clear workspace
clear j k z temp

fprintf('\nOptimization with discretization mesh conlcuded\n')

% clock stopped
fprintf('\nExecution concluded\n')
toc

%% PLOTTING

% if the user asks for plotting 
if options.plot

    fprintf('\n\nPlotting started\n')

    % plotting function
    run plotting.m
    
    fprintf('Plotting concluded\n')

end



% ASSIGNMENT 1
%   minimize the delta velocity required to depart from Jupiter and arrive
%   on Venus through a fly-by around the Earth

% TO DO
% - comment plotting.m
% - implement manual mesh for each window, apply multiple ga and the
% fmincon
% - comment main.m

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

% run configuration file, containing all numerical data
run config.m

%%  MESH OPTIMIZATION

% declare flyby modelling function as handle
JEV = @(dates) JEV_lamb_fb_lamb(dates, astro, options);

%
departure.window = linspace(departure.date_min, departure.date_max, optim.nmesh_dep);

for j = 1:1:size(departure.window,2);

    fly.window = linspace(departure.window(j) + departure.tof_min, ...
        departure.window(j) + departure.tof_max, optim.nmesh_fly);
    j

    for k = 1:1:size(fly.window,2);

        arrival.window = linspace(fly.window(k) + arrival.tof_min, ...
            fly.window(k) + arrival.tof_max, optim.nmesh_arr);

        for z = 1:1:size(arrival.window,2);

            temp.dates = [departure.window(j), fly.window(k), arrival.window(z)];
            temp.Dv = JEV(temp.dates);

            global results optim

            if temp.Dv < results.Dv_min
                results.dates = temp.dates;
                results.Dv_min = temp.Dv;
            end

            clear temp

        end
    end
end

clear j k z

fprintf('\nOptimization with discretization mesh conlcuded\n')

%% GA and FMINCON OPTIMIZATION

for index = 1:1:optim.noptim

    index
    
    % B constraint matrix for GA 
    optim.ga_Bcon = [ - results.dates(1) - departure.SW; ...
        - results.dates(2) - fly.SW; ...
        - results.dates(3) - arrival.SW; ...
        results.dates(1) + departure.SW; ...
        results.dates(2) + fly.SW; ...
        results.dates(3) + arrival.SW ];
    
    % optimization of Dv using ga
    [temp.dates, temp.Dv] = ...
        ga(JEV, optim.ga_nvars, optim.ga_Acon, optim.ga_Bcon, ...
        [], [], optim.ga_lb, [], [], options.ga_options);

    global results optim
    
%     % update fmincon guess according to ga result
%     fmincon_guess = temp.dates;
%     
%     % optimization of Dv using fmincon
%     [temp.dates, temp.Dv] = ...
%         fmincon(JEV, fmincon_guess, [], [], [], [], [], [], [], options.fmincon_options);

    if ~optim.ga_nanflag
        if temp.Dv < results.Dv_min
            results.dates = temp.dates;
            results.Dv_min = temp.Dv
        end
    else index = index - 1;
    end
    
    clear temp

end

 fprintf('\nGenetic algorithm and gradient method optimization conlcuded\n')

clear index

%-------------------------------------------------------------------------%


% clock stopped
fprintf('\nExecution concluded\n')
toc

%% PLOTTING


if options.plot
    fprintf('\nPlotting started\n')
    run plotting.m
    fprintf('\Plotting concluded\n')
end



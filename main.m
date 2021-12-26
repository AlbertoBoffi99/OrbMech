% ASSIGNMENT 1
%   minimize the delta velocity required to depart from Jupiter and arrive
%   on Venus through a fly-by around the Earth

% TO DO
% - implement manual mesh for each window, apply multiple ga on each window and the
% fmincon
% - see if plotting.m is optimized in terms of variables and datas
% - change plotting.m variables names
% - clear the workspace
% - comment every function

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

global out

% run configuration file, containing all numerical data
run config.m

%%  MESH OPTIMIZATION

% declare flyby modelling function as handle
JEV = @(dates) JEV_lamb_fb_lamb(dates, astro, options);

%
departure.window = linspace(departure.date_min, departure.date_max, optim.nmesh_dep);

temp.wb = waitbar(0, 'Pattern Optimization Running ...');

for j = 1:1:size(departure.window,2);

    fly.window = linspace(departure.window(j) + departure.tof_min, ...
        departure.window(j) + departure.tof_max, optim.nmesh_fly);
    
    waitbar(j/size(departure.window,2), temp.wb);

    for k = 1:1:size(fly.window,2);

        arrival.window = linspace(fly.window(k) + arrival.tof_min, ...
            fly.window(k) + arrival.tof_max, optim.nmesh_arr);

        for z = 1:1:size(arrival.window,2);

            temp.dates = [departure.window(j), fly.window(k), arrival.window(z)];
            temp.Dv = JEV(temp.dates);

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

clear j k z temp

fprintf('\nOptimization with discretization mesh conlcuded\n')

save('meshResults.mat', "results");

%% GA and FMINCON OPTIMIZATION

index = 1;
out.rp_norm = results.rp_norm;
 % non linear function for GA
out.ga_nonlincon = @(x) rp_nonlinFcn(x, astro.RE, astro.h_atm);
temp.wb = waitbar(0, 'GA Optimization Running ...');

for index = 1:1:optim.noptim

    waitbar(index/optim.noptim, temp.wb);
    
    % B constraint matrix for GA 
    optim.ga_Bcon = [ - results.dates(1) + departure.SW; ...
        - results.dates(2) + fly.SW; ...
        - results.dates(3) + arrival.SW; ...
        results.dates(1) + departure.SW; ...
        results.dates(2) + fly.SW; ...
        results.dates(3) + arrival.SW ];
    
    % optimization of Dv using GA
    [temp.dates, temp.Dv] = ...
        ga(JEV, optim.ga_nvars, optim.ga_Acon, optim.ga_Bcon, ...
        [], [], [], [], out.ga_nonlincon, [], options.ga_options);

    out.nanflag
    results
    temp

%     optim.fmincon_nonlincon = @(x) rp_nonlinFcn(x, out.rp_norm, astro.RE, astro.h_atm);
    
%     % update fmincon guess according to ga result
%     fmincon_guess = temp.dates;
%     
%     % optimization of Dv using fmincon
%     [temp.dates, temp.Dv] = ...
%         fmincon(JEV, fmincon_guess, optim.ga_Acon, optim.ga_Bcon, ...
%         [], [], [], [], optim.fmincon_nonlincon, options.fmincon_options);

    if temp.Dv < results.Dv_min
        results.dates = temp.dates;
        results.Dv_min = temp.Dv;
        results.ViT1 = out.ViT1;
        results.ViT2 = out.ViT2;
        results.DvT1 = out.DvT1;
        results.DvT2 = out.DvT2;
        results.Dvfb = out.Dvfb;
        results.rp = out.rp;
        results.rp_norm = out.rp_norm;
        results.vp_p = out.vp_p;
        results.vp_m = out.vp_m;
    end

end

 fprintf('\nGenetic algorithm and gradient method optimization conlcuded\n')

clear index temp

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



% FBTIME FILE
% Evaluates the time spent inside Earth's SOI

%-------------------------------------------------------------------------%

%% INPUT

% auxiliary multipurpose non-dimnesional variable
temp.n = 1000;
% time-span to plot the fly-by figure 
temp.plot_fbtime = 100000;
% contour plot cutting planes
temp.nplanes = 5;

%% SOI TIME

tspanH = linspace(0, temp.plot_fbtime, temp.n);
% outgoing hyperpola propagation
[timeH1, stateH1] = ode113( @(t,s) tbp_ode(t, s, astro.muE), tspanH, [results.rp,results.vp_p], options.ode_options);
%incoming hyperbola propagation
[timeH2, stateH2] = ode113( @(t,s) tbp_ode(t, s, astro.muE), -tspanH, [results.rp,results.vp_m], options.ode_options);
% extracting radii vectors
r_H1 = stateH1(:,1:3);
r_H2 = stateH2(:,1:3);
% SOI radius
r_soi = astroConstants(2)*((astroConstants(13)/astroConstants(1))/(astroConstants(4)/astroConstants(1)))^(2/5);

% ToF of 1st and 2nd hyperbolas
kk = 1;
while norm(r_H1(kk,:)) < r_soi && kk < length(tspanH)
    kk = kk + 1;
end
tofh_in = timeH1(kk);
kk = 1;
while norm(r_H2(kk,:)) < r_soi && kk < length(tspanH)
    kk = kk + 1;
end
tofh_out = -timeH2(kk);

% time spent inside the Earth's SOI
results.Dtfb = tofh_in + tofh_out;
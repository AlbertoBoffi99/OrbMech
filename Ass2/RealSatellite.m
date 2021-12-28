%% ATLAS 5 CENTAUR R/B NORAD 28473

kep0_real = [Kep_real(1,1), Kep_real(1,2), deg2rad(Kep_real(1,3)), deg2rad(Kep_real(1,4)), deg2rad(Kep_real(1,5)), 0];

% time interval of integration 
tspan = linspace(MJD2000*24*3600, MJD2000*24*3600 + 3*365*3600*24, 40000);

% integration solver options
options_sol = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);  

% propagation of Keplerian elements
[ time_Kep_real_prop, Kep_real_prop ] = ode113( @(t,s) gauss_propagation(t, s, muE, Re, J2, cr, A2m), tspan, kep0_real', options_sol );

%%

k_par = T*T/(3*365*3600*24)*length(tspan);
Kep_filt = movmean(Kep_real_prop(:,1), k_par);

figure()
% plot ((tspan-MJD2000*24*3600)/T, Kep(:,1));
% hold on
plot ((tspan-MJD2000*24*3600)/T, (Kep_real_prop(:,1)), (tspan_real-MJD2000*24*3600)/T, Kep_real(:,1), (tspan-MJD2000*24*3600)/T, Kep_filt);
grid on
legend ('GaussProp', 'Real')
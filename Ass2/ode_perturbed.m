function ds = ode_perturbed (t, s, muP, Re, J2, cr, A2m)

rr = s(1:3);
vv = s(4:6);

x= rr(1); y=rr(2); z=rr(3);
r = norm(rr);

%  J2 acceleration in cartesian (xyz)
aJ2_xyz = [  (3/2*J2*muP*Re^2/r^4)*(x/r*(5*z^2/r^2-1));
           (3/2*J2*muP*Re^2/r^4)*(y/r*(5*z^2/r^2-1));
           (3/2*J2*muP*Re^2/r^4)*(z/r*(5*z^2/r^2-3))];

% SRP acceleration in cartesian
AU = astroConstants(2);
pSR = 4.5e-6; % [N/m^2] SRP @ 1AU
[kep, muS] = uplanet(t/(24*60*60), 3);
[rE, vE] = kep2car( kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), muS);
r_sc2s = -(rE + rr);
aSRP_xyz = (- pSR * AU^2 * cr * A2m / (norm(r_sc2s))^3 * r_sc2s)/1e+3;

% gravitational acceleration in cartesian
a_xyz = [-muP/r^3*rr(1); -muP/r^3*rr(2); -muP/r^3*rr(3)];

% total acceleration
a_tot_xyz = a_xyz + aJ2_xyz + aSRP_xyz;

ds = [ vv(1); vv(2); vv(3); a_tot_xyz];

function dy = gauss_propagation(t,y,muP,Rp,J2,cr,A2m)

% Parameter's instant value
a = y(1); e = y(2); i = y(3); OM = y(4); om = y(5); th = y(6);
[r,v] = kep2car(a,e,i,OM,om,th,muP);

% Tangential-normal-outofplane RF
tt = v/norm(v);
hh = cross(r,v)/norm(cross(r,v));
nn = cross(hh,tt);

R_car2tnh = [tt'; nn'; hh']; % Rotation matrix from cartesian to TNH

% J2 related acceleration
r_norm = norm(r);
a_J2_car = ((3*J2*muP*Rp^2)/(2*r_norm^4))* [ (r(1)/r_norm)*(5*(r(3)^2/r_norm^2) - 1); 
                                             (r(2)/r_norm)*(5*(r(3)^2/r_norm^2) - 1); 
                                             (r(3)/r_norm)*(5*(r(3)^2/r_norm^2) - 3)];
% SRP
psr = 4.5e-6; % [N/m^2]
AU = astroConstants(2);
[kepE, muS] = uplanet(t/(24*3600), 3);
[rE, vE] = kep2car( kepE(1), kepE(2), kepE(3), kepE(4), kepE(5), kepE(6), muS);
r_sc_S = - (rE + r);
a_srp_car =  (- psr * (AU)^2 * cr * A2m / norm(r_sc_S)^3 * r_sc_S)/1e+3;

a_car = a_J2_car + a_srp_car;


a_tnh = R_car2tnh * a_car; % Acceleration in TNH frame
a_t = a_tnh(1);
a_n = a_tnh(2);
a_h = a_tnh(3);

% Parameters for Gauss planetary equations
b = a * sqrt(1-e^2);
%p = b^2/a;
n = sqrt(muP/a^3);
h = n*a*b;
th_star = th + om;

%Gauss planetary equations
r_n = norm(r);
v_n = norm(v);
a_dot  = (2*(a^2)*v_n/muP) * a_t;
e_dot  = (1/v_n) * (2*(e+cos(th))*a_t - (r_n/a)*sin(th)*a_n);
i_dot  = (r_n*cos(th_star)/h) * a_h;
OM_dot = (r_n*sin(th_star)/(h*sin(i))) * a_h;
om_dot = (1/(e*v_n)) * ( 2*sin(th)*a_t + (2*e + (r_n/a)*cos(th))*a_n ) - (r_n*sin(th_star)*cos(i)/(h*sin(i)))*a_h;
th_dot = h/r_n^2 - (1/(e*v_n)) * ( 2*sin(th)*a_t + (2*e + (r_n/a)*cos(th))*a_n );

% Derivative of the state
dy = [a_dot; e_dot; i_dot; OM_dot; om_dot; th_dot];

end
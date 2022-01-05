function ds = gauss_disturb(t, s, accfun, mu);
  
    %% INPUT
    a = s(1);
    e = s(2);
    i = s(3);
    RAAN = s(4);
    om = s(5);
    f = s(6);

    %% ODE
    
    b = a*sqrt(1-e^2);
    p = b^2/a;
    r_norm = p/(1+e*cos(f));
    n = sqrt(mu/a^3);
    v = sqrt(2*mu/r_norm - mu/a);
    h = n*a*b;
    fs = f + om;

    [r, ~] = kep2car(a, e, i, RAAN, om, f, mu);

    acc = accfun(r,t);

    % 3 trivial rotation matrixes
    Rot_RAAN = [cos(RAAN) sin(RAAN) 0; ...
              -sin(RAAN) cos(RAAN) 0; ...
              0 0 1];
    
    Rot_i = [1 0 0; ...
             0 cos(i) sin(i); ...
             0 -sin(i) cos(i)];
    
    Rot_omega = [cos(om) sin(om) 0; ...
              -sin(om) cos(om) 0; ...
              0 0 1];
    % rotation matrix from Geocentric to peri-focal 
    Rot_GE2PF = Rot_omega * Rot_i * Rot_RAAN;
    % radial velocity
    vr = mu/h*e*sin(f);
    % tangential velocity
    vtheta = mu/h*(1+e*cos(f));
    % flight path angle
    gamma = atan2(vr, vtheta);
    % rotation matrix from perifocal to NTH frame
    Rot_PF2NTH = [cos(gamma) sin(gamma) 0; ...
              -sin(gamma) cos(gamma) 0; ...
              0 0 1];
    % rotation matrix from Geocentric to NTH frame
    Rot_GE2NTH = Rot_PF2NTH * Rot_GE2PF;

    acc = Rot_GE2NTH*acc;

    an = acc(1);
    at = acc(2);
    ah = acc(3);

    ds = [ 2*a^2*v/mu * at;
        1/v*(2*(e+cos(f))*at - r_norm/a*sin(f)*an);
        r_norm*cos(fs)/h * ah;
        r_norm*sin(fs)/h/sin(i) * ah;
        1/e/v*(2*sin(f)*at + (2*e+r_norm/a*cos(f))*an)-(r_norm*sin(fs)*cos(i))/(h*sin(i)) * ah;
        h/r_norm^2 - 1/e/v*(2*sin(f)*at + (2*e+r_norm/a*cos(f))*an)];

end
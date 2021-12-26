function [c, ceq] = rp_nonlinFcn(x, RE, h_atm)

    global out

    if out.rp_norm < RE + h_atm
        c = abs(out.rp_norm - (RE + h_atm));
    else
        c = -1;
    end

    ceq = 0;

end
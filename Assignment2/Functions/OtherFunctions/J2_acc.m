function p = J2_acc(r);
  
  muP = astroConstants(13);
  Re = astroConstants(23);
  J2 = astroConstants(9);
  r_norm = sqrt(dot(r,r));

    p = [(3/2*J2*muP*Re^2/r_norm^4)*(r(1)/r_norm*(5*r(3)^2/r_norm^2-1));
       (3/2*J2*muP*Re^2/r_norm^4)*(r(2)/r_norm*(5*r(3)^2/r_norm^2-1));
       (3/2*J2*muP*Re^2/r_norm^4)*(r(3)/r_norm*(5*r(3)^2/r_norm^2-3))];

end

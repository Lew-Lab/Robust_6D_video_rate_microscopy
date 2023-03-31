function M2 = angle2M_2(theta,phi,gamma)

mux = cos(phi).*sin(theta);
muy = sin(phi).*sin(theta);
muz = cos(theta);

muxx = gamma.*mux.^2+(1-gamma)/3;
muyy = gamma.*muy.^2+(1-gamma)/3;
muzz = gamma.*muz.^2+(1-gamma)/3;
muxy = gamma.*mux.*muy;
muxz = gamma.*mux.*muz;
muyz = gamma.*muz.*muy;

M2 = [muxx;muyy;muzz;muxy;muxz;muyz];

end
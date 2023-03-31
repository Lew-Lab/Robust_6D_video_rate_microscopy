function [mux,muy,muz] = angle2xyz(theta,phi)

mux = cos(phi).*sin(theta);
muy = sin(phi).*sin(theta);
muz = cos(theta);

end
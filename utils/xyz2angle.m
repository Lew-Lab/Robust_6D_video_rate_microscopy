function [theta,phi] = xyz2angle(x,y,z)
theta = acos(z);
phi = atan2(y,x);

end
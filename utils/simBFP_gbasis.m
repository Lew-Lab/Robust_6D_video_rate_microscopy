function [gxPol,gyPol] = simBFP_gbasis(Microscopy)
lambda = Microscopy.wavelength; n1 = Microscopy.n1; NA = Microscopy.NA;  n2 = Microscopy.n2;
N = Microscopy.sampling_size; z = Microscopy.z; z2 = Microscopy.z2;

dx_true = Microscopy.pix_size/Microscopy.Magnitude*10^-9;
dx = n1*dx_true;
[eta,xi] = meshgrid(linspace(-1/(2*dx),1/(2*dx),N),linspace(-1/(2*dx),1/(2*dx),N));

xBFP = lambda*eta;
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho_max = NA/n1;

A = 1/sqrt((pi/3)*(4+(rho_max^2-4)*sqrt(1-rho_max^2)));

gxx = sqrt(1./sqrt(1-rho.^2)).*(cos(phi).^2.*sqrt(1-rho.^2)+sin(phi).^2)*A; 
gyy = gxx';
gxy = sqrt(1./sqrt(1-rho.^2)).*(sin(phi).*cos(phi).*(sqrt(1-rho.^2)-1))*A;
gyx = gxy';
gxz = sqrt(1./sqrt(1-rho.^2)).*rho.*cos(phi)*A;
gyz = gxz';

gxx(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gxy(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gxz(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gyx(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gyy(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gyz(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;

gxPol = cat(3,gxx,gxy,gxz);
gyPol = cat(3,gyx,gyy,gyz);

end

function f = circ(r)
f = (1 - sign(r-1))/2;
end
function [Ex_bfp,Ey_bfp] = simDipole(Microscopy)

lambda = Microscopy.wavelength; n1 = Microscopy.n1; NA = Microscopy.NA;  n2 = Microscopy.n2;
N = Microscopy.sampling_size; z = Microscopy.z; z2 = Microscopy.z2;

dx_true = Microscopy.pix_size/Microscopy.Magnitude*10^-9;
dx = n1*dx_true;
[eta,xi] = meshgrid(linspace(-1/(2*dx),1/(2*dx),N),linspace(-1/(2*dx),1/(2*dx),N));

xBFP = lambda*eta;
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho_max = NA/n1;

k1 = n1*(2*pi/lambda);
k2 = n2*(2*pi/lambda);
theta1 = asin(rho);
theta2 = asin((n1/n2)*sin(theta1));

%Fresnel coefficients
tp = 2*n2*cos(theta2)./(n2*cos(theta1) + n1*cos(theta2));
ts = 2*n2*cos(theta2)./(n1*cos(theta1) + n2*cos(theta2));

Esx = ts.*(cos(theta1)./cos(theta2)).*(n1/n2).*(-sin(phi));
Esy = ts.*(cos(theta1)./cos(theta2)).*(n1/n2).*cos(phi);
Epx = tp.*(n1/n2).*cos(phi).*cos(theta1);
Epy = tp.*(n1/n2).*sin(phi).*cos(theta1);
Epz = tp.*(-sin(theta1).*(n1/n2)^2.*(cos(theta1)./cos(theta2)));

% Exx - Ex contributed by mux
Exx = (1./sqrt(cos(theta1))).*(cos(phi).*Epx - sin(phi).*Esx).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
Exy = (1./sqrt(cos(theta1))).*(cos(phi).*Epy - sin(phi).*Esy).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
Exz = (1./sqrt(cos(theta1))).*(cos(phi).*Epz).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
Eyx = (1./sqrt(cos(theta1))).*(cos(phi).*Esx + sin(phi).*Epx).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
Eyy = (1./sqrt(cos(theta1))).*(cos(phi).*Esy + sin(phi).*Epy).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
Eyz = (1./sqrt(cos(theta1))).*(sin(phi).*Epz).*exp(1i*k1*z*cos(theta1));

Exx(rho >= rho_max) = 0;
Exy(rho >= rho_max) = 0;
Exz(rho >= rho_max) = 0;
Eyx(rho >= rho_max) = 0;
Eyy(rho >= rho_max) = 0;
Eyz(rho >= rho_max) = 0;

Ex_bfp = cat(3,Exx,Exy,Exz);
Ey_bfp = cat(3,Eyx,Eyy,Eyz);

end
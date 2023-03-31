%simulation of a dipole at an interface--By Adam Backer, Jan 1 2013
%Adapted from Axelrod, Journal of Microscopy article, 2012. (see journal club slides for more details)
function [basisBFPx_in,basisBFPy_in,ExBFP_in,EyBFP_in,basisBFPx_out,basisBFPy_out,ExBFP_out,EyBFP_out] = simDipole_BFP_v6(z,z2,zh,pmask,N,wavelength,n1,n2,nh,NA,M,pix_size,xy_ind)

% Yiyang modified based on 
% Tingting's simDipole_v5 and Oumeng's simDipole_v6

% simDipole_BFP_v1 is to generate the electric field at BFP and for further
% calculation and applying the Vortex Waveplate (raPol)

%inputs:
% x, y, z -- molecule coordinates note: z defocus length (unit:m)
% z2 -- distance of emitter below interface
% zh -- thickness of a thin film
% pmask -- phase mask
% N -- FT size
% n1 -- imaging media r.i.
% nh -- thin film r.i. (match this to n2, and set zh to zero if you don't want a thin film)
% n2 -- sample r.i. 
% N-- size
% NA -- numerical aperture
% M -- magnification

%outputs:
% basisBFPx_in -- x channel BFP images before phase mask modulation
% basisBFPy_in -- y channel BFP images before phase mask modulation
% ExBFP_in -- x channel electric field at BFP before phase mask modulation
% EyBFP_in -- y channel electric field at BFP before phase mask modulation
% basisBFPx_out -- x channel BFP images after phase mask modulation
% basisBFPy_out -- y channel BFP images after phase mask modulation
% ExBFP_out -- x channel electric field at BFP after phase mask modulation
% EyBFP_out -- y channel electric field at BFP after phase mask modulation

%simulation parameters%
lambda = wavelength;%wavelength

% parameter required to be added;


%calculate both pupil and image plane sampling, 
%one will affect the other, so make sure not to introduce aliasing

dx_true = (pix_size*1e-9/M);%image plane sampling
dx = n1*dx_true;%due to Abbe sine condition, scale by imaging medium r.i. (see appendix of my journal club)
%%???????  



dv = 1/(N*dx);%pupil sampling, related to image plane by FFT
% recall comb function, 1/dx is the preriod in fourier space(pupil space)
%%%????

%define pupil coordinates
temp=linspace((-1/(2*dx)),(1/(2*dx)),N);
[eta,xi] = meshgrid(temp);
% [eta,xi] = meshgrid(((-1/(2*dx))+(1/(2*N*dx))):dv:(-(1/(2*N*dx))+(1/(2*dx))),...
%     ((-1/(2*dx))+(1/(2*N*dx))):dv:(-(1/(N*2*dx))+(1/(2*dx))));

xBFP = lambda*eta;  % why scaled by wavelength??
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho_max = NA/n1;%pupil region of support determined by NA and imaging medium r.i.

k1 = n1*(2*pi/lambda);
kh = nh*(2*pi/lambda);
k2 = n2*(2*pi/lambda);
%rho(rho >= rho_max) = 0;
theta1 = asin(rho);%theta in matched medium
thetah = asin((n1/nh)*sin(theta1));%theta in thin film
theta2 = asin((n1/n2)*sin(theta1));%theta in mismatched medium
theta2 = real(theta2)-1i*abs(imag(theta2));

%%%%%%%%% Start
%Fresnel coefficients
tp_2h = 2*n2*cos(theta2)./(n2*cos(thetah) + nh*cos(theta2));
ts_2h = 2*n2*cos(theta2)./(nh*cos(thetah) + n2*cos(theta2));
tp_h1 = 2*nh*cos(thetah)./(nh*cos(theta1) + n1*cos(thetah));
ts_h1 = 2*nh*cos(thetah)./(n1*cos(theta1) + nh*cos(thetah));

rp_2h = (n2*cos(theta2) - nh*cos(thetah))./(n2*cos(theta2)+ nh*cos(thetah));
rs_2h = (nh*cos(theta2) - n2*cos(thetah))./(nh*cos(theta2)+ n2*cos(thetah));
rp_h1 = (nh*cos(thetah) - n1*cos(theta1))./(nh*cos(thetah)+ n1*cos(theta1));
rs_h1 = (n1*cos(thetah) - nh*cos(theta1))./(n1*cos(thetah)+ nh*cos(theta1));

%Axelrod's equations for E-fields in back focal plane

tp = tp_2h.*tp_h1.*exp(1i*kh*cos(thetah)*zh)./(1 + rp_2h.*rp_h1.*exp(2i*kh*zh*cos(thetah))); 
ts = ts_2h.*ts_h1.*exp(1i*kh*cos(thetah)*zh)./(1 + rs_2h.*rs_h1.*exp(2i*kh*zh*cos(thetah)));

% Es = ts.*(cos(theta1)./cos(theta2)).*(n1/n2).*(muy.*cos(phi) - mux.*sin(phi));
% Ep = tp.*((n1/n2).*(mux.*cos(phi) + muy.*sin(phi)).*cos(theta1) - muz*sin(theta1).*(n1/n2)^2.*(cos(theta1)./cos(theta2)));
% Esx - Es contributed by mux

% still has problem on how to get these to equation????????

%based on the equation above, seperating the mux,muy,muz compoment from two
%polirized electric field
Esx = ts.*(cos(theta1)./cos(theta2)).*(n1/n2).*(-sin(phi));
Esy = ts.*(cos(theta1)./cos(theta2)).*(n1/n2).*cos(phi);
Epx = tp.*(n1/n2).*cos(phi).*cos(theta1);
Epy = tp.*(n1/n2).*sin(phi).*cos(theta1);
Epz = tp.*(-sin(theta1).*(n1/n2)^2.*(cos(theta1)./cos(theta2)));

% Exx - Ex contributed by mux
% the first x represents x channel and y channel on the camera, the second
% x,y,z represents the compoment of mux, muy, muz from the orientation of
% dipole
Exx = (1./sqrt(cos(theta1))).*(cos(phi).*Epx - sin(phi).*Esx).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*z2*cos(theta2)); %added defocus aberration + depth aberration
Exy = (1./sqrt(cos(theta1))).*(cos(phi).*Epy - sin(phi).*Esy).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*z2*cos(theta2)); %added defocus aberration + depth aberration
Exz = (1./sqrt(cos(theta1))).*(cos(phi).*Epz).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*z2*cos(theta2)); %added defocus aberration + depth aberration
Eyx = (1./sqrt(cos(theta1))).*(cos(phi).*Esx + sin(phi).*Epx).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*z2*cos(theta2));
Eyy = (1./sqrt(cos(theta1))).*(cos(phi).*Esy + sin(phi).*Epy).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*z2*cos(theta2));
Eyz = (1./sqrt(cos(theta1))).*(sin(phi).*Epz).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*z2*cos(theta2));

% remove the electric component that is outside the accecptant region of
% objective lens
Exx(rho >= rho_max) = 0;
Exy(rho >= rho_max) = 0;
Exz(rho >= rho_max) = 0;
Eyx(rho >= rho_max) = 0;
Eyy(rho >= rho_max) = 0;
Eyz(rho >= rho_max) = 0;

ExBFP_in = cat(3,Exx,Exy,Exz);
EyBFP_in = cat(3,Eyx,Eyy,Eyz);

%% Oumeng's coordinate flipping
if xy_ind~=1
Exx = rot90(Exx);
Exy = rot90(Exy);
Exz = rot90(Exz);
Eyx = fliplr(Eyx);
Eyy = fliplr(Eyy);
Eyz = fliplr(Eyz);
end

%% Apply the phase mask to generate the electric field after phase modulating
[~,~,sizePmask] = size(pmask);
if sizePmask==2
    pmaskx=pmask(:,:,1);
    pmasky=pmask(:,:,2);
elseif sizePmask==1
    pmaskx=pmask;
    pmasky=pmask;
end

if xy_ind==1
    pmaskx=rot90(pmask,3);
    pmasky=pmask;
end
%for propagation from BFP E-field to image plane via tube-lens, paraxial
%approximation is in force.

Exx_out = Exx.*pmaskx;
Eyx_out = Eyx.*pmasky;
Exy_out = Exy.*pmaskx;
Eyy_out = Eyy.*pmasky;
Exz_out = Exz.*pmaskx;
Eyz_out = Eyz.*pmasky;
ExBFP_out = cat(3,Exx_out,Exy_out,Exz_out);
EyBFP_out = cat(3,Eyx_out,Eyy_out,Eyz_out);


basisBFPx_in(:,:,1) = abs(Exx).^2;
basisBFPx_in(:,:,2) = abs(Exy).^2;
basisBFPx_in(:,:,3) = abs(Exz).^2;
basisBFPx_in(:,:,4) = 2*real(Exx.*conj(Exy));
basisBFPx_in(:,:,5) = 2*real(Exx.*conj(Exz));
basisBFPx_in(:,:,6) = 2*real(Exy.*conj(Exz));
basisBFPy_in(:,:,1) = abs(Eyx).^2;
basisBFPy_in(:,:,2) = abs(Eyy).^2;
basisBFPy_in(:,:,3) = abs(Eyz).^2;
basisBFPy_in(:,:,4) = 2*real(Eyx.*conj(Eyy));
basisBFPy_in(:,:,5) = 2*real(Eyx.*conj(Eyz));
basisBFPy_in(:,:,6) = 2*real(Eyy.*conj(Eyz));

basisBFPx_out(:,:,1) = abs(Exx_out).^2;
basisBFPx_out(:,:,2) = abs(Exy_out).^2;
basisBFPx_out(:,:,3) = abs(Exz_out).^2;
basisBFPx_out(:,:,4) = 2*real(Exx_out.*conj(Exy_out));
basisBFPx_out(:,:,5) = 2*real(Exx_out.*conj(Exz_out));
basisBFPx_out(:,:,6) = 2*real(Exy_out.*conj(Exz_out));
basisBFPy_out(:,:,1) = abs(Eyx_out).^2;
basisBFPy_out(:,:,2) = abs(Eyy_out).^2;
basisBFPy_out(:,:,3) = abs(Eyz_out).^2;
basisBFPy_out(:,:,4) = 2*real(Eyx_out.*conj(Eyy_out));
basisBFPy_out(:,:,5) = 2*real(Eyx_out.*conj(Eyz_out));
basisBFPy_out(:,:,6) = 2*real(Eyy_out.*conj(Eyz_out));


end

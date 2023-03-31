function [basisImagex,basisImagey,basisBFPx,basisBFPy,basisImagex_dx,basisImagey_dx,basisImagex_dy,basisImagey_dy] = simDipole_raPol_v1(z,z2,zh,pmask,N,wavelength,n1,n2,nh,NA,M,pix_size,xy_ind)

% Yiyang modified based on simDipole_v5

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
% basisImagex -- x channel basis images
% basisImagey -- y channel basis images

%%simulation parameters--feel free to change%%
lambda = wavelength;%wavelength

%calculate both pupil and image plane sampling, 
%one will affect the other, so make sure not to introduce aliasing

dx_true = (pix_size*1e-9/M);%image plane sampling
dx = n1*dx_true;%due to Abbe sine condition, scale by imaging medium r.i. (see appendix of my journal club) 



dv = 1/(N*dx);%pupil sampling, related to image plane by FFT
% recall comb function, 1/dx is the preriod in fourier space(pupil space)
%%%????

%define pupil coordinates
temp=linspace((-1/(2*dx)),(1/(2*dx)),N);
[eta,xi] = meshgrid(temp);
% you can choose the radius of BFP arbitrarily for 2158C front system, 80 is used for 2158B back SLM imaging system

xBFP = lambda*eta;
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho_max = NA/n1;


%%
theta1 = asin(rho);
theta2 = asin((n1/n2)*sin(theta1));
a = cos(theta2)./cos(theta1);
%Fresnel coefficients
ts = 2*n1./(n1+n2*a);
tpxy = 2*n1./(n1+n2./a);
tpz = 2*n1./(n2+n1*a);

Esx = -sin(phi).*ts;
Esy = cos(phi).*ts;
Epx = cos(phi).*cos(theta1).*tpxy;
Epy = sin(phi).*cos(theta1).*tpxy;
Epz = -sin(theta1).*(n1/n2).*tpz;

% Exx - Ex contributed by mux
Exx = (1./sqrt(cos(theta1))).*(cos(phi).*Epx - sin(phi).*Esx);
Exy = (1./sqrt(cos(theta1))).*(cos(phi).*Epy - sin(phi).*Esy);
Exz = (1./sqrt(cos(theta1))).*(cos(phi).*Epz);
Eyx = (1./sqrt(cos(theta1))).*(cos(phi).*Esx + sin(phi).*Epx);
Eyy = (1./sqrt(cos(theta1))).*(cos(phi).*Esy + sin(phi).*Epy);
Eyz = (1./sqrt(cos(theta1))).*(sin(phi).*Epz);

Exx(rho >= rho_max) = 0;
Exy(rho >= rho_max) = 0;
Exz(rho >= rho_max) = 0;
Eyx(rho >= rho_max) = 0;
Eyy(rho >= rho_max) = 0;
Eyz(rho >= rho_max) = 0;

Ex_bfp = cat(3,Exx,Exy,Exz);
Ey_bfp = cat(3,Eyx,Eyy,Eyz);

%% Jones Vector to rotate electric field
%for propagation from BFP E-field to image plane via tube-lens, paraxial
%approximation is in force.
% CHIDO Par = pi
chidoPar = pi;
J11tmp = cos(chidoPar*n1/NA*rho/2) + 1j*sin(chidoPar*n1/NA*rho/2).*cos(phi);
J12tmp = -1j*sin(chidoPar*n1/NA*rho/2).*sin(phi);
J21tmp = -1j*sin(chidoPar*n1/NA*rho/2).*sin(phi);
J22tmp = cos(chidoPar*n1/NA*rho/2) - 1j*sin(chidoPar*n1/NA*rho/2).*cos(phi);
J11 = J11tmp*(1+1j) + J21tmp*(1-1j);
J12 = J12tmp*(1+1j) + J22tmp*(1-1j);
J21 = J11tmp*(1-1j) + J21tmp*(1+1j);
J22 = J12tmp*(1-1j) + J22tmp*(1+1j);

for i = 1:3
    Ex_bfp_out(:,:,i) = Ex_bfp(:,:,i).*J11+Ey_bfp(:,:,i).*J12;
    Ey_bfp_out(:,:,i) = Ex_bfp(:,:,i).*J21+Ey_bfp(:,:,i).*J22;
end

bfpExx = Ex_bfp_out(:,:,1); bfpExy = Ex_bfp_out(:,:,2); bfpExz = Ex_bfp_out(:,:,3);
bfpEyx = Ey_bfp_out(:,:,1); bfpEyy = Ey_bfp_out(:,:,2); bfpEyz = Ey_bfp_out(:,:,3);

%% Transform to image plane

imgExx = fftshift(fft2(bfpExx));
imgEyx = fftshift(fft2(bfpEyx));
imgExy = fftshift(fft2(bfpExy));
imgEyy = fftshift(fft2(bfpEyy));
imgExz = fftshift(fft2(bfpExz));
imgEyz = fftshift(fft2(bfpEyz));

% euqation from backer's paper Eq.22
basisImagex(:,:,1) = abs(imgExx).^2;
basisImagex(:,:,2) = abs(imgExy).^2;
basisImagex(:,:,3) = abs(imgExz).^2;
basisImagex(:,:,4) = 2*real(conj(imgExx).*imgExy);
basisImagex(:,:,5) = 2*real(conj(imgExx).*imgExz);
basisImagex(:,:,6) = 2*real(conj(imgExy).*imgExz);  
%the results are same with the above equation

basisImagey(:,:,1) = abs(imgEyx).^2;
basisImagey(:,:,2) = abs(imgEyy).^2;
basisImagey(:,:,3) = abs(imgEyz).^2;
basisImagey(:,:,4) = 2*real(imgEyx.*conj(imgEyy));
basisImagey(:,:,5) = 2*real(imgEyx.*conj(imgEyz));
basisImagey(:,:,6) = 2*real(imgEyy.*conj(imgEyz));

basisBFPx(:,:,1) = abs(bfpExx).^2;
basisBFPx(:,:,2) = abs(bfpExy).^2;
basisBFPx(:,:,3) = abs(bfpExz).^2;
basisBFPx(:,:,4) = 2*real(bfpExx.*conj(bfpExy));
basisBFPx(:,:,5) = 2*real(bfpExx.*conj(bfpExz));
basisBFPx(:,:,6) = 2*real(bfpExy.*conj(bfpExz));
basisBFPy(:,:,1) = abs(bfpEyx).^2;
basisBFPy(:,:,2) = abs(bfpEyy).^2;
basisBFPy(:,:,3) = abs(bfpEyz).^2;
basisBFPy(:,:,4) = 2*real(bfpEyx.*conj(bfpEyy));
basisBFPy(:,:,5) = 2*real(bfpEyx.*conj(bfpEyz));
basisBFPy(:,:,6) = 2*real(bfpEyy.*conj(bfpEyz));


% Gradient towards x and y direction
if xy_ind==1
    [X,~] = meshgrid(1*1e-9*(1:N),1e-9*(1:N));
    maskdx = exp(1j*2*pi*X/(N*dx_true));
    maskdy = rot90(maskdx);
    %maskdy = (maskdx);
    
    imgExx_dx = fftshift(fft2(bfpExx.*maskdx));
    imgEyx_dx = fftshift(fft2(bfpEyx.*maskdx));
    imgExy_dx = fftshift(fft2(bfpExy.*maskdx));
    imgEyy_dx = fftshift(fft2(bfpEyy.*maskdx));
    imgExz_dx = fftshift(fft2(bfpExz.*maskdx));
    imgEyz_dx = fftshift(fft2(bfpEyz.*maskdx));
    
    imgExx_dy = fftshift(fft2(bfpExx.*maskdy));
    imgEyx_dy = fftshift(fft2(bfpEyx.*maskdy));
    imgExy_dy = fftshift(fft2(bfpExy.*maskdy));
    imgEyy_dy = fftshift(fft2(bfpEyy.*maskdy));
    imgExz_dy = fftshift(fft2(bfpExz.*maskdy));
    imgEyz_dy = fftshift(fft2(bfpEyz.*maskdy));
    
    
    basisImagex_dx(:,:,1) = abs(imgExx_dx).^2;
    basisImagex_dx(:,:,2) = abs(imgExy_dx).^2;
    basisImagex_dx(:,:,3) = abs(imgExz_dx).^2;
    basisImagex_dx(:,:,4) = 2*real(conj(imgExx_dx).*imgExy_dx);
    basisImagex_dx(:,:,5) = 2*real(conj(imgExx_dx).*imgExz_dx);
    basisImagex_dx(:,:,6) = 2*real(conj(imgExy_dx).*imgExz_dx);  

    basisImagey_dx(:,:,1) = abs(imgEyx_dx).^2;
    basisImagey_dx(:,:,2) = abs(imgEyy_dx).^2;
    basisImagey_dx(:,:,3) = abs(imgEyz_dx).^2;
    basisImagey_dx(:,:,4) = 2*real(imgEyx_dx.*conj(imgEyy_dx));
    basisImagey_dx(:,:,5) = 2*real(imgEyx_dx.*conj(imgEyz_dx));
    basisImagey_dx(:,:,6) = 2*real(imgEyy_dx.*conj(imgEyz_dx));
    
    basisImagex_dy(:,:,1) = abs(imgExx_dy).^2;
    basisImagex_dy(:,:,2) = abs(imgExy_dy).^2;
    basisImagex_dy(:,:,3) = abs(imgExz_dy).^2;
    basisImagex_dy(:,:,4) = 2*real(conj(imgExx_dy).*imgExy_dy);
    basisImagex_dy(:,:,5) = 2*real(conj(imgExx_dy).*imgExz_dy);
    basisImagex_dy(:,:,6) = 2*real(conj(imgExy_dy).*imgExz_dy);  


    basisImagey_dy(:,:,1) = abs(imgEyx_dy).^2;
    basisImagey_dy(:,:,2) = abs(imgEyy_dy).^2;
    basisImagey_dy(:,:,3) = abs(imgEyz_dy).^2;
    basisImagey_dy(:,:,4) = 2*real(imgEyx_dy.*conj(imgEyy_dy));
    basisImagey_dy(:,:,5) = 2*real(imgEyx_dy.*conj(imgEyz_dy));
    basisImagey_dy(:,:,6) = 2*real(imgEyy_dy.*conj(imgEyz_dy));
 
else
    basisImagex_dx=[];
    basisImagey_dx=[];
    basisImagex_dy=[];
    basisImagey_dy=[];
end

end
function [basisImagex,basisImagey,basisBFPx,basisBFPy,basisImagex_dx,basisImagey_dx,basisImagex_dy,basisImagey_dy] = simDipole_raMVR_v1(z,zm,zh,N,wavelength,n1,n2,nh,NA,M,pixSize,bfpRadius,xy_ind,raRatio)

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

%Jputs:
% basisImagex -- x channel basis images
% basisImagey -- y channel basis images

%%simulation parameters--feel free to change%%
lambda = wavelength;%wavelength

%calculate both pupil and image plane sampling, 
%one will affect the other, so make sure not to introduce aliasing

dx_true = (pixSize*1e-9/M);%image plane sampling
dx = n1*dx_true;%due to Abbe sine condition, scale by imaging medium r.i. (see appendix of my journal club)

dv = 1/(N*dx);%pupil sampling, related to image plane by FFT
% recall comb function, 1/dx is the preriod in fourier space(pupil space)

%define pupil coordinates
temp=linspace((-1/(2*dx)),(1/(2*dx)),N);
[eta,xi] = meshgrid(temp);
% you can choose the radius of BFP arbitrarily for 2158C front system, 80 is used for 2158B back SLM imaging system

xBFP = lambda*eta;
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho_max = NA/n1;

%% Electric field calculation
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
Exx = (1./sqrt(cos(theta1))).*(cos(phi).*Epx - sin(phi).*Esx).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*zm*cos(theta2)); %added defocus aberration + depth aberration
Exy = (1./sqrt(cos(theta1))).*(cos(phi).*Epy - sin(phi).*Esy).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*zm*cos(theta2)); %added defocus aberration + depth aberration
Exz = (1./sqrt(cos(theta1))).*(cos(phi).*Epz).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*zm*cos(theta2)); %added defocus aberration + depth aberration
Eyx = (1./sqrt(cos(theta1))).*(cos(phi).*Esx + sin(phi).*Epx).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*zm*cos(theta2));
Eyy = (1./sqrt(cos(theta1))).*(cos(phi).*Esy + sin(phi).*Epy).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*zm*cos(theta2));
Eyz = (1./sqrt(cos(theta1))).*(sin(phi).*Epz).*exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*zm*cos(theta2));


Exx(rho >= rho_max) = 0;
Exy(rho >= rho_max) = 0;
Exz(rho >= rho_max) = 0;
Eyx(rho >= rho_max) = 0;
Eyy(rho >= rho_max) = 0;
Eyz(rho >= rho_max) = 0;

Ex_bfp = cat(3,Exx,Exy,Exz);
Ey_bfp = cat(3,Eyx,Eyy,Eyz);

%% raMVR modulation
%for propagation from BFP E-field to image plane via tube-lens, paraxial
%approximation is in force

bfpN = ceil(bfpRadius*2/rho_max);
if rem(bfpN,2) == 0
    bfpN = bfpN + 1;
end
[x,y] = meshgrid(linspace(-1,1,bfpN),linspace(-1,1,bfpN));
[t,r] = cart2pol(x,y);

phi_zeroPad = zeros(N);
phi_zeroPad((N-bfpN)/2+(1:bfpN),(N-bfpN)/2+(1:bfpN)) = t;

% Oumeng's calculation, different from mine
% Jone's Matrix to rotate electric field
J11 = cos(phi_zeroPad);
J22 = cos(phi_zeroPad);
J12 = sin(phi_zeroPad);
J21 = -sin(phi_zeroPad);
xChPart = zeros(N);
xChPart((N-bfpN)/2+(1:bfpN),(N-bfpN)/2+(1:bfpN)) = (x>0).*(y>0);
yChPart = zeros(N);
yChPart((N-bfpN)/2+(1:bfpN),(N-bfpN)/2+(1:bfpN)) = (t>pi/4).*(t<pi*3/4);

% J11 = cos(phi);
% J22 = cos(phi);
% J12 = sin(phi);
% J21 = -sin(phi);
% xChPart = (xBFP>0).*(yBFP>0);
% yChPart = (phi>pi/4).*(phi<pi*3/4);

for i = 1:3
    Ex_bfp_J(:,:,i) = Ex_bfp(:,:,i).*J11+Ey_bfp(:,:,i).*J12;
    Ey_bfp_J(:,:,i) = Ex_bfp(:,:,i).*J21+Ey_bfp(:,:,i).*J22;
end


for i = 1:3
    Ex_bfp_out(:,:,i) = [Ex_bfp_J(:,:,i).*xChPart,Ex_bfp_J(:,:,i).*rot90(xChPart),...
        Ex_bfp_J(:,:,i).*rot90(xChPart,2),Ex_bfp_J(:,:,i).*rot90(xChPart,3)].*sqrt(raRatio);
    Ey_bfp_out(:,:,i) = [Ey_bfp_J(:,:,i).*yChPart,Ey_bfp_J(:,:,i).*rot90(yChPart),...
        Ey_bfp_J(:,:,i).*rot90(yChPart,2),Ey_bfp_J(:,:,i).*rot90(yChPart,3)];
end

bfpExx = Ex_bfp_out(:,:,1); bfpExy = Ex_bfp_out(:,:,2); bfpExz = Ex_bfp_out(:,:,3);
bfpEyx = Ey_bfp_out(:,:,1); bfpEyy = Ey_bfp_out(:,:,2); bfpEyz = Ey_bfp_out(:,:,3);


%% Transform to image plane

for i = 1:3
    Ex_out(:,:,i) = [fftshift(fft2(Ex_bfp_J(:,:,i).*xChPart)),fftshift(fft2(Ex_bfp_J(:,:,i).*rot90(xChPart))),...
        fftshift(fft2(Ex_bfp_J(:,:,i).*rot90(xChPart,2))),fftshift(fft2(Ex_bfp_J(:,:,i).*rot90(xChPart,3)))]*sqrt(raRatio);
    Ey_out(:,:,i) = [fftshift(fft2(Ey_bfp_J(:,:,i).*yChPart)),fftshift(fft2(Ey_bfp_J(:,:,i).*rot90(yChPart))),...
        fftshift(fft2(Ey_bfp_J(:,:,i).*rot90(yChPart,2))),fftshift(fft2(Ey_bfp_J(:,:,i).*rot90(yChPart,3)))];
end

imgExx = Ex_out(:,:,1); imgExy = Ex_out(:,:,2); imgExz = Ex_out(:,:,3);
imgEyx = Ey_out(:,:,1); imgEyy = Ey_out(:,:,2); imgEyz = Ey_out(:,:,3);

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
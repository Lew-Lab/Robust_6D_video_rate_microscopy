function [B,E,B_dx,B_dy] = computeBasis(Microscopy,PSF)


lambda = Microscopy.wavelength; n1 = Microscopy.n1; NA = Microscopy.NA;  n = Microscopy.n2;
N = Microscopy.sampling_size; z = Microscopy.z;
zeroPadSize = N;
dx_true = Microscopy.pix_size/Microscopy.Magnitude*10^-9;
dx = n1*dx_true;
[eta,xi] = meshgrid(linspace(-1/(2*dx),1/(2*dx),N),linspace(-1/(2*dx),1/(2*dx),N));

x = lambda*eta;
y = lambda*xi;
[t,r] = cart2pol(x,y);    
    


t_zeroPad= t;
r_zeroPad = r;
x_zeroPad = x;
y_zeroPad = y;


switch PSF
    case 'xyPol' % xyPol
        J11 = ones(zeroPadSize);
        J22 = ones(zeroPadSize);
        J12 = zeros(zeroPadSize);
        J21 = zeros(zeroPadSize);
       [B,E,B_dx,B_dy] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);

    case 'raPol' % raPol
        J11 = cos(t_zeroPad);
        J22 = cos(t_zeroPad);
        J12 = sin(t_zeroPad);
        J21 = -sin(t_zeroPad);
        [B,E,B_dx,B_dy] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);

    case 'CHIDO_pi' % CHIDO Par = pi
        chidoPar = pi;
        J11tmp = cos(chidoPar*n/NA*r_zeroPad/2) + 1j*sin(chidoPar*n/NA*r_zeroPad/2).*cos(t_zeroPad);
        J12tmp = -1j*sin(chidoPar*n/NA*r_zeroPad/2).*sin(t_zeroPad);
        J21tmp = -1j*sin(chidoPar*n/NA*r_zeroPad/2).*sin(t_zeroPad);
        J22tmp = cos(chidoPar*n/NA*r_zeroPad/2) - 1j*sin(chidoPar*n/NA*r_zeroPad/2).*cos(t_zeroPad);
        J11 = J11tmp*(1+1j) + J21tmp*(1-1j);
        J12 = J12tmp*(1+1j) + J22tmp*(1-1j);
        J21 = J11tmp*(1-1j) + J21tmp*(1+1j);
        J22 = J12tmp*(1-1j) + J22tmp*(1+1j);
        [B,E,B_dx,B_dy] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);
    case 'CHIDO_halfpi' % CHIDO Par = pi/2
        chidoPar = pi/2;
        J11tmp = cos(chidoPar*n/NA*r_zeroPad/2) + 1j*sin(chidoPar*n/NA*r_zeroPad/2).*cos(t_zeroPad);
        J12tmp = -1j*sin(chidoPar*n/NA*r_zeroPad/2).*sin(t_zeroPad);
        J21tmp = -1j*sin(chidoPar*n/NA*r_zeroPad/2).*sin(t_zeroPad);
        J22tmp = cos(chidoPar*n/NA*r_zeroPad/2) - 1j*sin(chidoPar*n/NA*r_zeroPad/2).*cos(t_zeroPad);
        J11 = J11tmp*(1+1j) + J21tmp*(1-1j);
        J12 = J12tmp*(1+1j) + J22tmp*(1-1j);
        J21 = J11tmp*(1-1j) + J21tmp*(1+1j);
        J22 = J12tmp*(1-1j) + J22tmp*(1+1j);
        [B,E,B_dx,B_dy] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);
    case 'TS' % TS
        J12 = zeros(zeroPadSize);
        J21 = zeros(zeroPadSize);
        J11 = nan(zeroPadSize);
        phase_ramp = 3*pi;
        for i = 1:zeroPadSize
            for j = 1:zeroPadSize
                if (abs(j-zeroPadSize/2-0.5)>0.3*bfpRadius && abs(i-j)<bfpRadius) || i+j<zeroPadSize+1-bfpRadius || i+j > zeroPadSize+1+bfpRadius
                    J11(i,j) = exp(1j*phase_ramp/bfpRadius*(-0.5*(i-zeroPadSize/2-0.5)+0.5*sqrt(3)*abs(j-zeroPadSize/2-0.5)));
                else
                    J11(i,j) = exp(1j*phase_ramp/bfpRadius*(i-zeroPadSize/2-0.5));
                end
            end
        end
        J22 = rot90(J11);
        [B,E] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);
    case 'DH' % DH
        J12 = zeros(zeroPadSize);
        J21 = zeros(zeroPadSize);
        J11 = zeros(zeroPadSize);
        J11(zeroPadSize/2+[-127:128],zeroPadSize/2+[-127:128]) = ...
            exp(1j*double(imread('DH_biasAdjusted_stan_rescaled-resized_PupilRadius_80_Rot_90_xyShift_0_0.bmp'))/255*pi*2);
        J22 = rot90(J11);
        [B,E] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);
    case 'BS' % bisected
        J12 = zeros(zeroPadSize);
        J21 = zeros(zeroPadSize);
        J11 = nan(zeroPadSize);
        phase_ramp = 3*pi;
        J11 = exp(1j*phase_ramp/bfpRadius*abs(x_zeroPad/2*bfpSampling));
        J22 = rot90(J11);
        [B,E] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);
    case 'Astigma' % astigmatic
        z_ast = 200e-9;
        J12 = zeros(zeroPadSize);
        J21 = zeros(zeroPadSize);
        J11 = exp(1j*2*pi*n/lambda*z_ast*cos(asin(abs(x_zeroPad))));
        J22 = rot90(J11);
        [B,E] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);
    case 'Biplane' % biplane
        z_bip = 350e-9;
        J12 = zeros(zeroPadSize);
        J21 = zeros(zeroPadSize);
        J11 = exp(1j*2*pi*n/lambda*z_bip/2*cos(asin(sqrt(x_zeroPad.^2+y_zeroPad.^2))));
        J22 = exp(-1j*2*pi*n/lambda*z_bip/2*cos(asin(sqrt(x_zeroPad.^2+y_zeroPad.^2))));
        [B,E] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);
    case 'Vortex'
        J12 = zeros(zeroPadSize);
        J21 = zeros(zeroPadSize);
        J11 = exp(1j*t_zeroPad);
        J22 = rot90(J11);
        [B,E,B_dx,B_dy] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);
        
    case 'optimized'
        J12 = zeros(zeroPadSize);
        J21 = zeros(zeroPadSize);
        J11 = zeros(zeroPadSize);
        J11(zeroPadSize/2+[-127:128],zeroPadSize/2+[-127:128]) = ...
            exp(1j*double(imread('polar_CRB_angle_v3_v2.bmp'))/255*pi*2);
        J11 = rot90(J11,Microscopy.rot);
        J22 = rot90(J11);
        [B,E] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22);

end

end




function [B,E,B_dx,B_dy] = computeBasis_BFP_mismatch(Microscopy,J11,J12,J21,J22)

NA = Microscopy.NA;
n = Microscopy.n1;
rm = NA/n;
n2 = Microscopy.n2;

bfpRadius = Microscopy.bfp_radius;
bfpSampling = ceil(bfpRadius*2/rm);
if rem(bfpSampling,2) == 1
    bfpSampling = bfpSampling + 1;
end
img_sz = Microscopy.image_size;

zeroPadSize = Microscopy.sampling_size;  %=lambda*f/d_pixel_size
if rem(zeroPadSize,2) == 1
    zeroPadSize = zeroPadSize + 1;
end

[Ex_bfp,Ey_bfp] = simDipole(Microscopy);
PSFregion = zeroPadSize/2+[-(img_sz-1)/2:(img_sz-1)/2];

for i = 1:3
    Ex_bfp_J(:,:,i) = Ex_bfp(:,:,i).*J11+Ey_bfp(:,:,i).*J12;
    Ey_bfp_J(:,:,i) = Ex_bfp(:,:,i).*J21+Ey_bfp(:,:,i).*J22;
    Ex(:,:,i) = fftshift(fft2(Ex_bfp_J(:,:,i)));
    Ey(:,:,i) = fftshift(fft2(Ey_bfp_J(:,:,i)));
        
end

E = [Ex(PSFregion,PSFregion,:),Ey(PSFregion,PSFregion,:)];

%E = [Ex,Ey];
Ex_tmp = Ex;
Ey_tmp = Ey;


for i = 1:3
    Bx_tmp(:,:,i) = abs(Ex_tmp(:,:,i)).^2;
end
Bx_tmp(:,:,4) = 2*real(Ex_tmp(:,:,1).*conj(Ex_tmp(:,:,2)));
Bx_tmp(:,:,5) = 2*real(Ex_tmp(:,:,1).*conj(Ex_tmp(:,:,3)));
Bx_tmp(:,:,6) = 2*real(Ex_tmp(:,:,2).*conj(Ex_tmp(:,:,3)));

for i = 1:3
    By_tmp(:,:,i) = abs(Ey_tmp(:,:,i)).^2;
end
By_tmp(:,:,4) = 2*real(Ey_tmp(:,:,1).*conj(Ey_tmp(:,:,2)));
By_tmp(:,:,5) = 2*real(Ey_tmp(:,:,1).*conj(Ey_tmp(:,:,3)));
By_tmp(:,:,6) = 2*real(Ey_tmp(:,:,2).*conj(Ey_tmp(:,:,3)));

B_tmp = [Bx_tmp,By_tmp];
intensity = sum(sum(B_tmp(:,:,1)+B_tmp(:,:,2)+B_tmp(:,:,3)))/3;
B_tmp = [Bx_tmp(PSFregion,PSFregion,:),By_tmp(PSFregion,PSFregion,:)];
%intensity = sum(sum(B_tmp(:,:,1)+B_tmp(:,:,2)+B_tmp(:,:,3)))/3;
B = B_tmp/intensity/0.9051;


% for dx,dy

 if Microscopy.xy_ind==1
       [X,~] = meshgrid(1e-9*(1:zeroPadSize),1e-9*(1:zeroPadSize));
       maskdx = exp(1j*2*pi*X/(zeroPadSize*58.5e-9));
       maskdy = rot90(maskdx);
       for i=1:3
           Ex_dx(:,:,i) = fftshift(fft2(Ex_bfp_J(:,:,i).*maskdx));
           Ey_dx(:,:,i) = fftshift(fft2(Ey_bfp_J(:,:,i).*maskdx)); 
           Ex_dy(:,:,i) = fftshift(fft2(Ex_bfp_J(:,:,i).*maskdy));
           Ey_dy(:,:,i) = fftshift(fft2(Ey_bfp_J(:,:,i).*maskdy)); 
       end
       
       Ex_tmp_dx = Ex_dx;
       Ey_tmp_dx = Ey_dx;
       Ex_tmp_dy = Ex_dy;
       Ey_tmp_dy = Ey_dy;
       
       %dx
       for i = 1:3
            Bx_tmp_dx(:,:,i) = abs(Ex_tmp_dx(:,:,i)).^2;
       end
       Bx_tmp_dx(:,:,4) = 2*real(Ex_tmp_dx(:,:,1).*conj(Ex_tmp_dx(:,:,2)));
       Bx_tmp_dx(:,:,5) = 2*real(Ex_tmp_dx(:,:,1).*conj(Ex_tmp_dx(:,:,3)));
       Bx_tmp_dx(:,:,6) = 2*real(Ex_tmp_dx(:,:,2).*conj(Ex_tmp_dx(:,:,3)));
       
       for i = 1:3
           By_tmp_dx(:,:,i) = abs(Ey_tmp_dx(:,:,i)).^2;
       end
       By_tmp_dx(:,:,4) = 2*real(Ey_tmp_dx(:,:,1).*conj(Ey_tmp_dx(:,:,2)));
       By_tmp_dx(:,:,5) = 2*real(Ey_tmp_dx(:,:,1).*conj(Ey_tmp_dx(:,:,3)));
       By_tmp_dx(:,:,6) = 2*real(Ey_tmp_dx(:,:,2).*conj(Ey_tmp_dx(:,:,3)));
       
       %dy
       for i = 1:3
            Bx_tmp_dy(:,:,i) = abs(Ex_tmp_dy(:,:,i)).^2;
       end
       Bx_tmp_dy(:,:,4) = 2*real(Ex_tmp_dy(:,:,1).*conj(Ex_tmp_dy(:,:,2)));
       Bx_tmp_dy(:,:,5) = 2*real(Ex_tmp_dy(:,:,1).*conj(Ex_tmp_dy(:,:,3)));
       Bx_tmp_dy(:,:,6) = 2*real(Ex_tmp_dy(:,:,2).*conj(Ex_tmp_dy(:,:,3)));
       
       for i = 1:3
           By_tmp_dy(:,:,i) = abs(Ey_tmp_dy(:,:,i)).^2;
       end
       By_tmp_dy(:,:,4) = 2*real(Ey_tmp_dy(:,:,1).*conj(Ey_tmp_dy(:,:,2)));
       By_tmp_dy(:,:,5) = 2*real(Ey_tmp_dy(:,:,1).*conj(Ey_tmp_dy(:,:,3)));
       By_tmp_dy(:,:,6) = 2*real(Ey_tmp_dy(:,:,2).*conj(Ey_tmp_dy(:,:,3)));
       
       B_tmp_dx = [Bx_tmp_dx,By_tmp_dx];
       intensity = sum(sum(B_tmp_dx(:,:,1)+B_tmp_dx(:,:,2)+B_tmp_dx(:,:,3)))/3;
       B_tmp_dx = [Bx_tmp_dx(PSFregion,PSFregion,:),By_tmp_dx(PSFregion,PSFregion,:)];
       B_dx = B_tmp_dx/intensity/0.9051;
       
       B_tmp_dy = [Bx_tmp_dy,By_tmp_dx];
       intensity = sum(sum(B_tmp_dy(:,:,1)+B_tmp_dy(:,:,2)+B_tmp_dy(:,:,3)))/3;
       B_tmp_dy = [Bx_tmp_dy(PSFregion,PSFregion,:),By_tmp_dy(PSFregion,PSFregion,:)];
       B_dy = B_tmp_dy/intensity/0.9051;
       
        
  else
        B_dx=[];
        B_dy=[];
  end

end


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





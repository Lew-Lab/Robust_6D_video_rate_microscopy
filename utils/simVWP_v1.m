function [basis_matrix,basisBFPx_out,basisBFPy_out,Ex_bfp_out,Ey_bfp_out] = simVWP_v1(Microscopy,Ex_bfp,Ey_bfp)

lambda = Microscopy.wavelength; n1 = Microscopy.n1; NA = Microscopy.NA;  n2 = Microscopy.n2;
N = Microscopy.sampling_size; z = Microscopy.z;
zeroPadSize = N;
dx_true = Microscopy.pix_size/Microscopy.Magnitude*10^-9;
dx = n1*dx_true;
[eta,xi] = meshgrid(linspace(-1/(2*dx),1/(2*dx),N),linspace(-1/(2*dx),1/(2*dx),N));

rm = NA/n1;

bfpRadius = Microscopy.bfp_radius;
bfpSampling = ceil(bfpRadius*2/rm);
img_sz = Microscopy.image_size;

if rem(bfpSampling,2) == 1
    bfpSampling = bfpSampling + 1;
end

zeroPadSize = Microscopy.sampling_size;  %=lambda*f/d_pixel_size
if rem(zeroPadSize,2) == 1
    zeroPadSize = zeroPadSize + 1;
end

x = lambda*eta;
y = lambda*xi;
[t,r] = cart2pol(x,y);    

t_zeroPad= t;
r_zeroPad = r;
x_zeroPad = x;
y_zeroPad = y;

PSFregion = zeroPadSize/2+[-(img_sz-1)/2:(img_sz-1)/2];

% Jones Vector for raPol
J11 = cos(t_zeroPad);
J22 = cos(t_zeroPad);
J12 = sin(t_zeroPad);
J21 = -sin(t_zeroPad);

Ex_bfp_out = zeros(size(Ex_bfp));
Ey_bfp_out = zeros(size(Ey_bfp));
for i = 1:3
    Ex_bfp_out(:,:,i) = Ex_bfp(:,:,i).*J11+Ey_bfp(:,:,i).*J12;
    Ey_bfp_out(:,:,i) = Ex_bfp(:,:,i).*J21+Ey_bfp(:,:,i).*J22;
    Ex(:,:,i) = fftshift(fft2(Ex_bfp_out(:,:,i)));
    Ey(:,:,i) = fftshift(fft2(Ey_bfp_out(:,:,i)));
        
end

E = [Ex(PSFregion,PSFregion,:),Ey(PSFregion,PSFregion,:)];

Exx_out = Ex_bfp_out(:,:,1);
Exy_out = Ex_bfp_out(:,:,2);
Exz_out = Ex_bfp_out(:,:,3);
Eyx_out = Ey_bfp_out(:,:,1);
Eyy_out = Ey_bfp_out(:,:,2);
Eyz_out = Ey_bfp_out(:,:,3);

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

% Electric field and basis images at image plane
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
B = B_tmp/intensity;

basisx = Bx_tmp(PSFregion,PSFregion,:)./intensity./0.8019;
basisy = By_tmp(PSFregion,PSFregion,:)./intensity./0.8019;

% Reshape the image into a 2D N*6 matrix

for i = 1:6  
    A = reshape(basisx(:,:,i),img_sz*img_sz,1);
    B = reshape(basisy(:,:,i),img_sz*img_sz,1);
    basis_matrix(:,i) = cat(1,A,B);
    
end

end
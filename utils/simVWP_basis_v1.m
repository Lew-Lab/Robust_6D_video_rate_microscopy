function [basis_matrix,Bgradx,Bgrady] = simVWP_basis_v1(Microscopy)

%% Imaging system parameters
name=Microscopy.mask;
rot=Microscopy.rot;
lambda=Microscopy.wavelength;
bfp_radius=Microscopy.bfp_radius;
n1=Microscopy.n1;
n2=Microscopy.n2;
nh=Microscopy.nh;
NA=Microscopy.NA;
Magnitude=Microscopy.Magnitude;
sampling_size=Microscopy.sampling_size;%=lambda*f/d_pixel_size
image_size=Microscopy.image_size;
pix_size=Microscopy.pix_size;
upsamping=Microscopy.upsampling;
pixelSizeUpsampling=Microscopy.pixelSizeUpsampling;
zf = Microscopy.z;
z2 = Microscopy.z2;
xy_ind = Microscopy.xy_ind;
dx_true = pix_size/Magnitude*10^-9;
dx = n1*dx_true; 
[eta,xi] = meshgrid(linspace(-1/(2*dx),1/(2*dx),sampling_size),linspace(-1/(2*dx),1/(2*dx),sampling_size));
xBFP = lambda*eta;
yBFP = lambda*xi;

rho_max = NA/n1;
bfpSampling = ceil(bfp_radius*2/rho_max);

if rem(bfpSampling,2) == 1
    bfpSampling = bfpSampling + 1;
end


if rem(sampling_size,2) == 1
    sampling_size = sampling_size + 1;
end

pmask = nan;

%% Electrical Field at BFP
[Ex_bfp,Ey_bfp] = simDipole(Microscopy);

[phi,rho] = cart2pol(xBFP,yBFP);    

% t_zeroPad= phi;
% r_zeroPad = rho;
% x_zeroPad = x;
% y_zeroPad = y;

% phase pattern representing defocus/lateral shift 
pdx = exp(1j*2*pi*n1/lambda*1e-9*xBFP); % phase pattern corresponding to molecule shifting 1 nm along x
pdy = exp(1j*2*pi*n1/lambda*1e-9*yBFP); % phase pattern corresponding to molecule shifting 1 nm along y
pdz = exp(1j*2*pi*n1/lambda*1e-9*cos(asin(rho))); % phase pattern corresponding to molecule shifting 1 nm along z
pz = exp(1j*2*pi*n1/lambda*zf*cos(asin(rho)));
p = cat(3,ones(sampling_size).*pz,pdx.*pz,pdy.*pz,pdz.*pz); % group them so code is cleaner later - ones(N) represents no phase pattern

PSFregion = sampling_size/2+(-(image_size-1)/2:(image_size-1)/2);

% Jones Vector for raPol
J11 = cos(phi);
J22 = cos(phi);
J12 = sin(phi);
J21 = -sin(phi);

Ex_bfp_out = zeros(size(Ex_bfp));
Ey_bfp_out = zeros(size(Ey_bfp));
for i = 1:3
    Ex_bfp_out(:,:,i) = Ex_bfp(:,:,i).*J11+Ey_bfp(:,:,i).*J12;
    Ey_bfp_out(:,:,i) = Ex_bfp(:,:,i).*J21+Ey_bfp(:,:,i).*J22;
    Ex(:,:,i) = fftshift(fft2(Ex_bfp_out(:,:,i)));
    Ey(:,:,i) = fftshift(fft2(Ey_bfp_out(:,:,i)));
        
end

%% Basis Images - raPol

% Electric field and basis images at image plane
Ex_tmp = Ex;
Ey_tmp = Ey;

for i = 1:3
    Bx_tmp(:,:,i) = abs(Ex_tmp(:,:,i)).^2;
end
Bx_tmp(:,:,4) = real(Ex_tmp(:,:,1).*conj(Ex_tmp(:,:,2))+conj(Ex_tmp(:,:,1)).*Ex_tmp(:,:,2));
Bx_tmp(:,:,5) = real(Ex_tmp(:,:,1).*conj(Ex_tmp(:,:,3))+conj(Ex_tmp(:,:,1)).*Ex_tmp(:,:,3));
Bx_tmp(:,:,6) = real(Ex_tmp(:,:,2).*conj(Ex_tmp(:,:,3))+conj(Ex_tmp(:,:,2)).*Ex_tmp(:,:,3));

for i = 1:3
    By_tmp(:,:,i) = abs(Ey_tmp(:,:,i)).^2;
end
By_tmp(:,:,4) = real(Ey_tmp(:,:,1).*conj(Ey_tmp(:,:,2))+conj(Ey_tmp(:,:,1)).*Ey_tmp(:,:,2));
By_tmp(:,:,5) = real(Ey_tmp(:,:,1).*conj(Ey_tmp(:,:,3))+conj(Ey_tmp(:,:,1)).*Ey_tmp(:,:,3));
By_tmp(:,:,6) = real(Ey_tmp(:,:,2).*conj(Ey_tmp(:,:,3))+conj(Ey_tmp(:,:,2)).*Ey_tmp(:,:,3));

B_tmp = [Bx_tmp,By_tmp];
intensity = sum(sum(B_tmp(:,:,1)+B_tmp(:,:,2)+B_tmp(:,:,3)))./3;
B_tmp = [Bx_tmp(PSFregion,PSFregion,:),By_tmp(PSFregion,PSFregion,:)];
%intensity = sum(sum(B_tmp(:,:,1)+B_tmp(:,:,2)+B_tmp(:,:,3)))/3;

basisx = Bx_tmp(PSFregion,PSFregion,:)./intensity;
basisy = By_tmp(PSFregion,PSFregion,:)./intensity;

% Reshape the image into a 2D N*6 matrix

for i = 1:6  
    A = reshape(basisx(:,:,i),image_size*image_size,1);
    B = reshape(basisy(:,:,i),image_size*image_size,1);
    basis_matrix(:,i) = cat(1,A,B);
    
end

%% FT to find image plane electric field

for j = 1:4 % j = 1:4 represent no phase, shift along x, shift along y, shift along z, respectively
    for i = 1:3 % 3 basis fields
        Ex_bfpTmp(:,:,i) = Ex_bfp_out(:,:,i).*p(:,:,j);
        Ey_bfpTmp(:,:,i) = Ey_bfp_out(:,:,i).*p(:,:,j);
        Ex(:,:,i) = fftshift(fft2(Ex_bfpTmp(:,:,i)));
        Ey(:,:,i) = fftshift(fft2(Ey_bfpTmp(:,:,i)));
    end
    Elist{j} = [Ex(PSFregion,PSFregion,:),Ey(PSFregion,PSFregion,:)];
end

%% image plane basis
for j = 1:4 % j = 1:4 represent no phase, shift along x, shift along y, shift along z, respectively
    E_tmp = Elist{j};
    for i = 1:3
        B_tmp(:,:,i) = abs(E_tmp(:,:,i)).^2;
    end
    B_tmp(:,:,4) = 2*real(E_tmp(:,:,1).*conj(E_tmp(:,:,2)));
    B_tmp(:,:,5) = 2*real(E_tmp(:,:,1).*conj(E_tmp(:,:,3)));
    B_tmp(:,:,6) = 2*real(E_tmp(:,:,2).*conj(E_tmp(:,:,3)));
    B_tmp = B_tmp/sum(sum(B_tmp(:,:,1))); % normalization to 1
    Blist{j} = B_tmp;
end
Bgradx = (Blist{2}-Blist{1})*1e9; % gradient along x
Bgrady = (Blist{3}-Blist{1})*1e9; % gradient along y

end
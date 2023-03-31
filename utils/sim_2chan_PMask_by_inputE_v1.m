function [basisImgs,basisBFPx_out,basisBFPy_out,Ex_bfp_out,Ey_bfp_out] = sim_2chan_PMask_by_inputE_v1(Microscopy,pmask,Ex_BFP_inMask,Ey_BFP_inMask)

name=Microscopy.mask;
rot=Microscopy.rot;
wavelength=Microscopy.wavelength;
bfp_radius=Microscopy.bfp_radius;
n1=Microscopy.n1;
n2=Microscopy.n2;
nh=Microscopy.nh;
NA=Microscopy.NA;
Magnitude=Microscopy.Magnitude;
sampling_size=Microscopy.sampling_size;
image_size=Microscopy.image_size;
pix_size=Microscopy.pix_size;
upsamping=Microscopy.upsampling;
pixelSizeUpsampling=Microscopy.pixelSizeUpsampling;
zf = Microscopy.z;
z2 = Microscopy.z2;
xy_ind = Microscopy.xy_ind;
%MaskResized=mask_resize(name,bfp_radius,rot,sampling_size);

img_sz = Microscopy.image_size;
zeroPadSize = Microscopy.sampling_size;  %=lambda*f/d_pixel_size
if rem(zeroPadSize,2) == 1
    zeroPadSize = zeroPadSize + 1;
end
PSFregion = zeroPadSize/2+(-(img_sz-1)/2:(img_sz-1)/2);

%% Get the electric field from the input
Exx_inMask = Ex_BFP_inMask(:,:,1);
Eyx_inMask = Ey_BFP_inMask(:,:,1);
Exy_inMask = Ex_BFP_inMask(:,:,2);
Eyy_inMask = Ey_BFP_inMask(:,:,2);
Exz_inMask = Ex_BFP_inMask(:,:,3);
Eyz_inMask = Ey_BFP_inMask(:,:,3);

%% Generate the phase mask
for pmask_channel = 1:2
    angle_temp = pmask(:,:,pmask_channel);
    angle_temp = rot90(angle_temp,rot);
    angle_1_sgl_ch = zeros(sampling_size, sampling_size)*127;
    angle_temp = im2double(angle_temp,'indexed');
    angleResize = imresize(angle_temp,upsamping);
    
    center = round((sampling_size+1)/2);
    [psf_radius,~] = size(angleResize);
    region = center+(-psf_radius/2:psf_radius/2-1);
    angle_1_sgl_ch(region,region)=angleResize(:,:);
    angle_1(:,:,pmask_channel) = angle_1_sgl_ch;
end

% angle_1 = angle_1-min(angle_1,[],'all');
if isfield(Microscopy, 'zernikeCoeff')
   angle_1 = pmask_corrected(angle_1,Microscopy);
end
%angle_1 = imrotate(angle_1,5);
pmask = exp(1i*angle_1);

%  for i = 1:sampling_size
%     for j = 1:sampling_size
%         mask(i,j) = exp(1j*5*pi/80*abs(j-sampling_size/2-0.5));
%     end
% end


% apply the phase mask to generate the electric field after phase modulating
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

Exx_out = Exx_inMask.*pmaskx;
Eyx_out = Eyx_inMask.*pmasky;
Exy_out = Exy_inMask.*pmaskx;
Eyy_out = Eyy_inMask.*pmasky;
Exz_out = Exz_inMask.*pmaskx;
Eyz_out = Eyz_inMask.*pmasky;
Ex_bfp_out = cat(3,Exx_out,Exy_out,Exz_out);
Ey_bfp_out = cat(3,Eyx_out,Eyy_out,Eyz_out);

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


%% Generate the basis images with E field
for i = 1:3
    Ex(:,:,i) = fftshift(fft2(Ex_bfp_out(:,:,i)));
    Ey(:,:,i) = fftshift(fft2(Ey_bfp_out(:,:,i)));
        
end

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

% basisx = Bx_tmp(PSFregion,PSFregion,:)./intensity./0.8019;
% basisy = By_tmp(PSFregion,PSFregion,:)./intensity./0.8019;

basisx = Bx_tmp(PSFregion,PSFregion,:)./intensity;
basisy = By_tmp(PSFregion,PSFregion,:)./intensity;

% Reshape the image into a 2D N*6 matrix

for i = 1:6  
    A = reshape(basisx(:,:,i),img_sz*img_sz,1);
    B = reshape(basisy(:,:,i),img_sz*img_sz,1);
    basis_matrix(:,i) = cat(1,A,B);
    
end

basisImgs = basis_matrix;

end
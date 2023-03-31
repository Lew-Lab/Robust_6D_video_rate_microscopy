function [BasisMatrix,BFPMatrix,BasisMatrix_dx,BasisMatrix_dy] = computeBasisWithModulation(Microscopy,pmask,bfpMod)
name=Microscopy.mask;
rot=Microscopy.rot;
wavelength=Microscopy.wavelength;
bfpRadius=Microscopy.bfp_radius;
n1=Microscopy.n1;
n2=Microscopy.n2;
nh=Microscopy.nh;
NA=Microscopy.NA;
Magnitude=Microscopy.Magnitude;
sampling_size=Microscopy.sampling_size;
imgSize=Microscopy.image_size;
pixSize=Microscopy.pix_size;
upsamping=Microscopy.upsampling;
pixSizeUpsampling=Microscopy.pixelSizeUpsampling;
zf = Microscopy.z;
z2 = Microscopy.z2;
xyInd = Microscopy.xy_ind;
%MaskResized=mask_resize(name,bfp_radius,rot,sampling_size);

%% 1. Generate the mask
if isnan(pmask)
    angle_temp = imread(name);
    angle_temp = rot90(angle_temp,rot);
    angle_1 = ones(sampling_size, sampling_size)*127;
    angle_temp = im2double(angle_temp,"indexed");
    angleResize = imresize(angle_temp,upsamping);
    
    center = round((sampling_size+1)/2);
    [psf_radius, ~] = size(angleResize);
    region = center+(-psf_radius/2:psf_radius/2-1);
    angle_1(region ,region)=angleResize(:,:);
    angle_1 = ((angle_1/255)*2*pi)-pi;
else
    angle_1 = zeros(sampling_size, sampling_size);
    angle_temp = pmask;
    angleResize = imresize(angle_temp,upsamping);
    center = round((sampling_size+1)/2);
    [psf_radius, ~] = size(angleResize);
    region = center+(-psf_radius/2:psf_radius/2-1);
    angle_1(region ,region)=angleResize(:,:);    
end

angle_1 = angle_1-min(angle_1,[],"all");

%% 2. Calculate Basis with different Unitary Operator
switch bfpMod
    case "xyPol"
        angle_1 = zeros(sampling_size,sampling_size);
        mask = exp(1i*angle_1);
        [basisImgx,basisImgy,BFP_img_x,BFP_img_y,basisImgx_dx,basisImgy_dx,basisImgx_dy,basisImgy_dy] = simDipole_v5(zf,z2,0,mask,sampling_size,wavelength,n1,n2,nh,NA,Magnitude,pixSize,xyInd);

    case "pmask"
        mask = exp(1i*angle_1);
        [basisImgx,basisImgy,BFP_img_x,BFP_img_y,basisImgx_dx,basisImgy_dx,basisImgx_dy,basisImgy_dy] = simDipole_v5(zf,z2,0,mask,sampling_size,wavelength,n1,n2,nh,NA,Magnitude,pixSize,xyInd);
   
    case "raPol"
        angle_1 = zeros(sampling_size,sampling_size);
        mask = exp(1i*angle_1);
        [basisImgx,basisImgy,BFP_img_x,BFP_img_y,basisImgx_dx,basisImgy_dx,basisImgx_dy,basisImgy_dy] = simDipole_raPol_v1(zf,z2,0,mask,sampling_size,wavelength,n1,n2,nh,NA,Magnitude,pixSize,xyInd);

    case "MVR"
        raRatio = 1;
        [basisImgx,basisImgy,BFP_img_x,BFP_img_y,basisImgx_dx,basisImgy_dx,basisImgx_dy,basisImgy_dy] = simDipole_raMVR_v1(zf,z2,0,sampling_size,wavelength,n1,n2,nh,NA,Magnitude,pixSize,bfpRadius,xyInd,raRatio);


end

%% 3. Normalize the Basis Images
if bfpMod ~= "MVR"
    intensity = 1/3*(sum(sum(basisImgx(:,:,1)))+sum(sum(basisImgy(:,:,1))))+1/3*(sum(sum(basisImgx(:,:,2)))+sum(sum(basisImgy(:,:,2))))+1/3*(sum(sum(basisImgx(:,:,3)))+sum(sum(basisImgy(:,:,3))));
    region = round(sampling_size/2)+[-(imgSize-1)/2:(imgSize-1)/2];
    basisx = basisImgx(region,region,:);
    basisy = basisImgy(region,region,:);
    
    [a,b]=size(BFP_img_x);
    a = a-1; b = b-1;
    
    BFP_x = BFP_img_x(round(a/2)-bfpRadius+1:round(a/2)+bfpRadius,round(a/2)-bfpRadius+1:round(a/2)+bfpRadius,:);
    BFP_y = BFP_img_y(round(a/2)-bfpRadius+1:round(a/2)+bfpRadius,round(a/2)-bfpRadius+1:round(a/2)+bfpRadius,:);
    
    
    %intensity = 1/3*(sum(sum(basisx(:,:,1)))+sum(sum(basisy(:,:,1))))+1/3*(sum(sum(basisx(:,:,2)))+sum(sum(basisy(:,:,2))))+1/3*(sum(sum(basisx(:,:,3)))+sum(sum(basisy(:,:,3))));
    basisx = basisx./intensity; %0.9051 for 71 %0.8250 for 35  %0.8109 for 31%normaliza basis images
    basisy = basisy./intensity;  %0.9051 for 71
    
    
    if xyInd ==1
        basisx_dx = basisImgx_dx(region,region,:);
        basisy_dx = basisImgy_dx(region,region,:);
        basisx_dx = basisx_dx./intensity; 
        basisy_dx = basisy_dx./intensity; 
        
        basisx_dy = basisImgx_dy(region,region,:);
        basisy_dy = basisImgy_dy(region,region,:);
        basisx_dy = basisx_dy./intensity; 
        basisy_dy = basisy_dy./intensity; 
        
        BasisMatrix_dx = cat(2,basisx_dx,basisy_dx);
        BasisMatrix_dy = cat(2,basisx_dy,basisy_dy);
    else
       BasisMatrix_dx=[];
       BasisMatrix_dy=[];
        
    end
    intensity_BFP = 1/3*(sum(sum(BFP_x(:,:,1)))+sum(sum(BFP_y(:,:,1))))+1/3*(sum(sum(BFP_x(:,:,2)))+sum(sum(BFP_y(:,:,2))))+1/3*(sum(sum(BFP_x(:,:,3)))+sum(sum(BFP_y(:,:,3))));
    BFP_x = BFP_x./intensity_BFP;
    BFP_y = BFP_y./intensity_BFP;
    
    BasisMatrix = cat(2,basisx,basisy);
    BFPMatrix   = cat(2,BFP_x,BFP_y);
else
    intensity = 1/3*(sum(sum(basisImgx(:,:,1)))+sum(sum(basisImgy(:,:,1))))+1/3*(sum(sum(basisImgx(:,:,2)))+sum(sum(basisImgy(:,:,2))))+1/3*(sum(sum(basisImgx(:,:,3)))+sum(sum(basisImgy(:,:,3))));
    region = round(sampling_size/2)+[-(imgSize-1)/2:(imgSize-1)/2];

    basisx = zeros(imgSize,imgSize*4,6); basisy = zeros(imgSize,imgSize*4,6);

    for i = 1:4
        basisx_curCh = basisImgx(:,(sampling_size*(i-1)+1):sampling_size*i,:);
        basisx(:,(imgSize*(i-1)+1):imgSize*i,:) = basisx_curCh(region,region,:);
        basisy_curCh = basisImgy(:,(sampling_size*(i-1)+1):sampling_size*i,:);
        basisy(:,(imgSize*(i-1)+1):imgSize*i,:) = basisy_curCh(region,region,:);
    end

    basisx = basisx./intensity; 
    basisy = basisy./intensity;

    % BFP normalization for MVR did not finish
    BasisMatrix = cat(2,basisx,basisy);
    BFPMatrix = 0; 
    BasisMatrix_dx=0;BasisMatrix_dy = 0;

end


end
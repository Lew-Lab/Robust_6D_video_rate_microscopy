close all;clear;clc
addpath(genpath('utils'));
addpath(genpath('phasemask'));

%% Imaging Parameters and Basis Images
signal=5000;
background=2;

Microscopy=struct();
Microscopy.wavelength = 593e-9; 
Microscopy.n1 = 1.518;
Microscopy.n2 = 1.334;
% Microscopy.n2 = 1.518;
Microscopy.nh = 1.518;
Microscopy.NA = 1.4;
Microscopy.pixelSizeUpsampling = 1;
Microscopy.upsampling = 1;
Microscopy.pix_size=6500/Microscopy.pixelSizeUpsampling;
Microscopy.bfp_radius = 80*Microscopy.upsampling;
Microscopy.Magnitude = 111.1111;
dx = (Microscopy.pix_size*10^-9)/Microscopy.Magnitude;
Microscopy.sampling_size = round(Microscopy.wavelength*Microscopy.bfp_radius/dx/Microscopy.NA)*Microscopy.pixelSizeUpsampling;
if rem(Microscopy.sampling_size,2)==0
    Microscopy.sampling_size=Microscopy.sampling_size+1;
end
Microscopy.image_size = 45*Microscopy.pixelSizeUpsampling; %the final image is a 45 by 45 image, odd number
Microscopy.z = -500*10^-9; %(unit:m) focal length (negative value, e.g. -500*10^-9, objective focused 500nm above the coverslip);
Microscopy.z2 = [-300e-09:dx:1000e-09]; %(unit:m) SM's z position (axial location of SM, positive)
Microscopy.zh = -0*10^-9;
Microscopy.xy_ind = 1;
Microscopy.rot=0;
imgSz = Microscopy.image_size;
% import opt PSF  
pmask = nan; %if you want to use .mat as pmask, then load .mat data to pmask
%coordinate are same as the microscopy coordinate
Microscopy.mask='pixOL_v12_com_vShift_M1_hShift_0.bmp'; 

%% Calculate the basis images
temp = hot(300); map1(:,1) = min(temp(:,3),1); map1(:,2) = min(temp(:,1)*0.50+temp(:,2)*0.50,1); map1(:,3) = min(temp(:,1)*1,1);

% [basisMatrix,BFP,Bgradx,Bgrady] = computeBasisWithModulation(Microscopy,pmask,'pmask'); % simulate PSF with different BFP modulation
[basisMatrices,mask_opt,BFP] = basis_matrix_pixel_based_v3_in(Microscopy,pmask); % Tingting's calculations
[a,numBasImgs] = size(basisMatrices{1});
numPixel = sqrt(a/2);

%% Reshape basis images to match inverse estimator format

for i = 1:length(basisMatrices)
    basisMatrix = basisMatrices{i};
    [a,~] = size(basisMatrix);
    numPixel = sqrt(a/2);
    % Mxx
    temp = reshape(basisMatrix(:,1),numPixel,2*numPixel);
    basisXXx=temp(:,1:numPixel);  basisXXy=temp(:,numPixel+1:2*numPixel);
    basisXX = [basisXXx basisXXy]; Bstruct.Bxx{i} = basisXX;
    % Myy
    temp = reshape(basisMatrix(:,2),numPixel,2*numPixel);
    basisYYx=temp(:,1:numPixel);  basisYYy=temp(:,numPixel+1:2*numPixel);
    basisYY = [basisYYx basisYYy]; Bstruct.Byy{i} = basisYY;
    % Mzz
    temp = reshape(basisMatrix(:,3),numPixel,2*numPixel);
    basisZZx=temp(:,1:numPixel);  basisZZy=temp(:,numPixel+1:2*numPixel);
    basisZZ = [basisZZx basisZZy]; Bstruct.Bzz{i} = basisZZ;
    % Mxy
    temp = reshape(basisMatrix(:,4),numPixel,2*numPixel);
    basisXYx=temp(:,1:numPixel);  basisXYy=temp(:,numPixel+1:2*numPixel);
    basisXY = [basisXYx basisXYy]; Bstruct.Bxy{i} = basisXY;
    % Mxz
    temp = reshape(basisMatrix(:,5),numPixel,2*numPixel);
    basisXZx=temp(:,1:numPixel);  basisXZy=temp(:,numPixel+1:2*numPixel);
    basisXZ = [basisXZx basisXZy]; Bstruct.Bxz{i} = basisXZ;
    % Myz
    temp = reshape(basisMatrix(:,6),numPixel,2*numPixel);
    basisYZx=temp(:,1:numPixel);  basisYZy=temp(:,numPixel+1:2*numPixel);
    basisYZ = [basisYZx basisYZy]; Bstruct.Byz{i} = basisYZ;
end

save('BpixOL.mat', 'Bstruct', 'Microscopy');
%{
M_iso = [1/3;1/3;1/3;0;0;0]; % Isotropic dipole

theta = pi/2; phi = 0; gamma = 1; % dipole orientation and wobble parameter

M = angle2M_2(theta,phi,gamma);

I = signal*basisMatrix*M_iso + background;
IxPol = I(1:numPixel^2); IxPol = reshape(IxPol,numPixel,numPixel);
IyPol = I(numPixel^2+1:2*numPixel^2); IyPol = reshape(IyPol,numPixel,numPixel);
%}

%% Plot the basis images
%{
[a,~] = size(basisMatrix);
numPixel = sqrt(a/2);
% Mxx
temp = reshape(basisMatrix(:,1),numPixel,2*numPixel);
basisXXx=temp(:,1:numPixel);  basisXXy=temp(:,numPixel+1:2*numPixel);
% Myy
temp = reshape(basisMatrix(:,2),numPixel,2*numPixel);
basisYYx=temp(:,1:numPixel);  basisYYy=temp(:,numPixel+1:2*numPixel);
% Mzz
temp = reshape(basisMatrix(:,3),numPixel,2*numPixel);
basisZZx=temp(:,1:numPixel);  basisZZy=temp(:,numPixel+1:2*numPixel);
% Mxy
temp = reshape(basisMatrix(:,4),numPixel,2*numPixel);
basisXYx=temp(:,1:numPixel);  basisXYy=temp(:,numPixel+1:2*numPixel);
% Mxz
temp = reshape(basisMatrix(:,5),numPixel,2*numPixel);
basisXZx=temp(:,1:numPixel);  basisXZy=temp(:,numPixel+1:2*numPixel);
% Myz
temp = reshape(basisMatrix(:,6),numPixel,2*numPixel);
basisYZx=temp(:,1:numPixel);  basisYZy=temp(:,numPixel+1:2*numPixel);

% Plot the Image Plane
figure('Position',[0 0 800 600]);
subplot(2,3,1);
basisXX = [basisXXx;basisXXy];
imagesc(basisXX);
axis image; colorbar; axis off;
title("B_{xx}");

subplot(2,3,2);
basisYY = [basisYYx;basisYYy];
imagesc(basisYY);
axis image; colorbar; axis off;
title("B_{yy}");

subplot(2,3,3);
basisZZ = [basisZZx;basisZZy];
imagesc(basisZZ);
axis image; colorbar; axis off;
title("B_{zz}");

subplot(2,3,4);
basisXY = [basisXYx;basisXYy];
imagesc(basisXY);
axis image; colorbar; axis off;
title("B_{xy}");

subplot(2,3,5);
basisXZ = [basisXZx;basisXZy];
imagesc(basisXZ);
axis image; colorbar; axis off;
title("B_{xz}");

subplot(2,3,6);
basisYZ = [basisYZx;basisYZy];
imagesc(basisYZ);
axis image; colorbar; axis off;
title("B_{yz}");

%% Plot the simulated SM images

% Ground Truth
figure();
pos1 = [0.3 0.51 0.4 0.4];
ax1 = axes('Position',pos1);
subplot('Position',pos1);
imagesc(IxPol);
text(2,3,'X-Pol','FontWeight','bold','Color','white')
axis image; colormap(ax1,'hot'); axis off; colorbar;
% title('Basis Image: X-pol');
pos2 = [0.3 0.09 0.4 0.4];
ax2 = axes('Position',pos2);
subplot('Position',pos2);
imagesc(IyPol); hold on;
line([38,38+3.4188],[40,40],'LineWidth',3,'Color','white');
text(2,3,'Y-Pol','FontWeight','bold','Color','white')
text(38+3.4188/2,42,'200 nm','Color','white','HorizontalAlignment','center','FontWeight','bold');
axis image; colormap(ax2,map1); axis off; colorbar;
% title('Basis Image: Y-pol');
%}

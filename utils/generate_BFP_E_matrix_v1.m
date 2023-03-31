function [BFP_x_in,BFP_y_in,ExBFP_in,EyBFP_in,BFP_x_out,BFP_y_out,ExBFP_out,EyBFP_out] = generate_BFP_E_matrix_v1(Microscopy)

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

% 1. generate the mask
angle_temp = imread(name);
angle_temp = rot90(angle_temp,rot);
angle_1 = ones(sampling_size, sampling_size)*127;
angle_temp = im2double(angle_temp,'indexed');
angleResize = imresize(angle_temp,upsamping);

center = round((sampling_size+1)/2);
[psf_radus,~] = size(angleResize);
region = center+(-psf_radus/2:psf_radus/2-1);
angle_1(region ,region+3 )=angleResize(:,:);
angle_1 = ((angle_1/255)*2*pi)-pi;

angle_1 = angle_1-min(angle_1,[],'all');
if isfield(Microscopy, 'zernikeCoeff')
   angle_1 = pmask_corrected(angle_1,Microscopy);
end
%angle_1 = imrotate(angle_1,5);
mask = exp(1i*angle_1);

%  for i = 1:sampling_size
%     for j = 1:sampling_size
%         mask(i,j) = exp(1j*5*pi/80*abs(j-sampling_size/2-0.5));
%     end
% end

%2.a generate the basis image
[BFPImg_x_in,BFPImg_y_in,ExBFP_in,EyBFP_in,BFPImg_x_out,BFPImg_y_out,ExBFP_out,EyBFP_out] = simDipole_BFP_v6(zf,z2,0,mask,sampling_size,wavelength,n1,n2,nh,NA,Magnitude,pix_size,xy_ind);

[a,b]=size(BFPImg_x_in);
a = a-1;
b = b-1;
BFP_x_in = BFPImg_x_in(round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,:);
BFP_y_in = BFPImg_y_in(round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,:);

[a,b]=size(BFPImg_x_out);
a = a-1;
b = b-1;
BFP_x_out = BFPImg_x_out(round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,:);
BFP_y_out = BFPImg_y_out(round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,:);



intensity1 = 1/3*(sum(sum(BFP_x_in(:,:,1)))+sum(sum(BFP_y_in(:,:,1))))+1/3*(sum(sum(BFP_x_in(:,:,2)))+sum(sum(BFP_y_in(:,:,2))))+1/3*(sum(sum(BFP_x_in(:,:,3)))+sum(sum(BFP_y_in(:,:,3))));
BFP_x_in = BFP_x_in/intensity1;
BFP_y_in = BFP_y_in/intensity1;

intensity2 = 1/3*(sum(sum(BFP_x_out(:,:,1)))+sum(sum(BFP_y_out(:,:,1))))+1/3*(sum(sum(BFP_x_out(:,:,2)))+sum(sum(BFP_y_out(:,:,2))))+1/3*(sum(sum(BFP_x_out(:,:,3)))+sum(sum(BFP_y_out(:,:,3))));
BFP_x_out = BFP_x_out/intensity2;
BFP_y_out = BFP_y_out/intensity2;

%2.b reshape the image to (75*75+75*75) by 6 matrix
[BFP_size,~,~]=size(BFP_x_in);
BFP_matrix_in = zeros(BFP_size*BFP_size+BFP_size*BFP_size,6);
BFP_matrix_out = zeros(BFP_size*BFP_size+BFP_size*BFP_size,6);

for i = 1:6  
    A1 = reshape(BFP_x_in(:,:,i),BFP_size*BFP_size,1);
    B1 = reshape(BFP_y_in(:,:,i),BFP_size*BFP_size,1);
    BFP_matrix_in(:,i) = cat(1,A1,B1);

    A2 = reshape(BFP_x_out(:,:,i),BFP_size*BFP_size,1);
    B2 = reshape(BFP_y_out(:,:,i),BFP_size*BFP_size,1);
    BFP_matrix_out(:,i) = cat(1,A2,B2);
end


end


%% Built-in functions
function pmask = pmask_corrected(pmaskOrg,Microscopy)


lambda=Microscopy.wavelength;
NA=Microscopy.NA;
M=Microscopy.Magnitude;
N=Microscopy.sampling_size;
pix_size=Microscopy.pix_size;
n1=Microscopy.n1;
zernikeCoeff = Microscopy.zernikeCoeff;

rm = NA/n1;
dx_true = pix_size*1e-9/M;%image plane sampling
dx = n1*dx_true;
temp=linspace((-1/(2*dx)),(1/(2*dx)),N);
[eta,xi] = meshgrid(temp);
xBFP = lambda*eta; 
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho = rho/rm;

%calculate n_c,m_c

zernMax = sqrt(2*(length(zernikeCoeff)+1)+1/4)-1/2-1;
nindx=[];mindx=[];
for i=1:zernMax                   
    nindx=[nindx i*ones(1,mod(i+1,2)+2*ceil(i/2))];
    mindx=[mindx -i:2:i];
end
mindx(nindx<=1)=[];
zernikeCoeff(nindx<=1,:)=[];
nindx(nindx<=1)=[];

Z= zernfun(nindx,mindx,rho,phi); 
pmaskZern = reshape(reshape(Z,N*N,[])*zernikeCoeff,N,N,[]);
thred(:,:,1)=rho;
thred(:,:,2)=rho;
pmaskZern(thred>1)=0;

pmask = pmaskOrg+pmaskZern;
            
end
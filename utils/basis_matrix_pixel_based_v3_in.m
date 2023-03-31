function [basis_matrices,mask,BFP_matrix,basis_matrix_dx,basis_matrix_dy] = basis_matrix_pixel_based_v3_in(Microscopy,pmask)
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

xy_ind = Microscopy.xy_ind;
%MaskResized=mask_resize(name,bfp_radius,rot,sampling_size);
basis_matrices = cell(1, length(Microscopy.z2));
for j = 1:length(Microscopy.z2)
    z2 = Microscopy.z2(j);
    %1. generate the mask
    if isnan(pmask)==1
        angle_temp = imread(name);
        angle_temp = rot90(angle_temp,rot);
        angle_1 = ones(sampling_size, sampling_size)*127;
        angle_temp = im2double(angle_temp,'indexed');
        angleResize = imresize(angle_temp,upsamping);

        center = round((sampling_size+1)/2);
        [psf_radius, ~] = size(angleResize);
        region = center+(-psf_radius/2:psf_radius/2-1);
        angle_1(region ,region+3)=angleResize(:,:);
        angle_1 = ((angle_1/255)*2*pi)-pi;
    else
        angle_1 = ones(sampling_size, sampling_size);
        angle_temp = pmask;
        angleResize = imresize(angle_temp,upsamping);
        center = round((sampling_size+1)/2);
        [psf_radius, ~] = size(angleResize);
        region = center+(-psf_radius/2:psf_radius/2-1);
        angle_1(region ,region+3)=angleResize(:,:);

    end

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
    [basisImagex,basisImagey,BFP_image_x,BFP_image_y,basisImagex_dx,basisImagey_dx,basisImagex_dy,basisImagey_dy] = simDipole_v5(zf,z2,0,mask,sampling_size,wavelength,n1,n2,nh,NA,Magnitude,pix_size,xy_ind);
    intensity = 1/3*(sum(sum(basisImagex(:,:,1)))+sum(sum(basisImagey(:,:,1))))+1/3*(sum(sum(basisImagex(:,:,2)))+sum(sum(basisImagey(:,:,2))))+1/3*(sum(sum(basisImagex(:,:,3)))+sum(sum(basisImagey(:,:,3))));
    region = round(sampling_size/2)+[-(image_size-1)/2:(image_size-1)/2];
    basisx = basisImagex(region,region,:);
    basisy = basisImagey(region,region,:);

    [a,b]=size(BFP_image_x);
    a = a-1;
    b = b-1;

    BFP_x = BFP_image_x(round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,:);
    BFP_y = BFP_image_y(round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,round(a/2)-bfp_radius+1:round(a/2)+bfp_radius,:);


    %intensity = 1/3*(sum(sum(basisx(:,:,1)))+sum(sum(basisy(:,:,1))))+1/3*(sum(sum(basisx(:,:,2)))+sum(sum(basisy(:,:,2))))+1/3*(sum(sum(basisx(:,:,3)))+sum(sum(basisy(:,:,3))));
    basisx = basisx./intensity/0.6025; %0.9051 for 71 %0.8250 for 35  %0.8109 for 31%normaliza basis images
    basisy = basisy./intensity/0.6025;  %0.9051 for 71


    if xy_ind ==1
        basisx_dx = basisImagex_dx(region,region,:);
        basisy_dx = basisImagey_dx(region,region,:);
        basisx_dx = basisx_dx./intensity/0.6025; 
        basisy_dx = basisy_dx./intensity/0.6025; 

        basisx_dy = basisImagex_dy(region,region,:);
        basisy_dy = basisImagey_dy(region,region,:);
        basisx_dy = basisx_dy./intensity/0.6025; 
        basisy_dy = basisy_dy./intensity/0.6025; 

        for i = 1:6  
        A = reshape(basisx_dx(:,:,i),image_size*image_size,1);
        B = reshape(basisy_dx(:,:,i),image_size*image_size,1);
        basis_matrix_dx(:,i) = cat(1,A,B);


        A = reshape(basisx_dy(:,:,i),image_size*image_size,1);
        B = reshape(basisy_dy(:,:,i),image_size*image_size,1);
        basis_matrix_dy(:,i) = cat(1,A,B);

        end
    else
       basis_matrix_dx=[];
       basis_matrix_dy=[];

    end
    intensity2 = 1/3*(sum(sum(BFP_x(:,:,1)))+sum(sum(BFP_y(:,:,1))))+1/3*(sum(sum(BFP_x(:,:,2)))+sum(sum(BFP_y(:,:,2))))+1/3*(sum(sum(BFP_x(:,:,3)))+sum(sum(BFP_y(:,:,3))));
    BFP_x = BFP_x/intensity2;
    BFP_y = BFP_y/intensity2;


    %2.b reshape the image to (75*75+75*75) by 6 matrix
    [BFP_size,~,~]=size(BFP_x);
    basis_matrix = zeros(image_size*image_size+image_size*image_size,6);

    BFP_matrix = zeros(BFP_size*BFP_size+BFP_size*BFP_size,6);

    for i = 1:6  
        A = reshape(basisx(:,:,i),image_size*image_size,1);
        B = reshape(basisy(:,:,i),image_size*image_size,1);
        basis_matrix(:,i) = cat(1,A,B);


        A1 = reshape(BFP_x(:,:,i),BFP_size*BFP_size,1);
        B1 = reshape(BFP_y(:,:,i),BFP_size*BFP_size,1);
        BFP_matrix(:,i) = cat(1,A1,B1);


    end
    basis_matrices{j} = basis_matrix;
end

end



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
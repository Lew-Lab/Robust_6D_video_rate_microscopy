function [basisgxPol,basisgyPol,basisGxPol,basisGyPol] = simIP_Gbasis_noNormalize(Microscopy,pmask,BFPMod)

lambda = Microscopy.wavelength; n1 = Microscopy.n1; NA = Microscopy.NA;  n2 = Microscopy.n2;
z = Microscopy.z; z2 = Microscopy.z2; 
sampling_size = Microscopy.sampling_size; image_size = Microscopy.image_size; upsampling = Microscopy.upsampling;
name = Microscopy.mask; rot=Microscopy.rot;

dx_true = Microscopy.pix_size/Microscopy.Magnitude*10^-9;
dx = n1*dx_true;
[eta,xi] = meshgrid(linspace(-1/(2*dx),1/(2*dx),sampling_size),linspace(-1/(2*dx),1/(2*dx),sampling_size));

xBFP = lambda*eta;
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho_max = NA/n1;

%1. generate the mask
if isnan(pmask)
    angle_temp = imread(name);
    angle_temp = rot90(angle_temp,rot);
    angle_1 = zeros(sampling_size, sampling_size);
    angle_temp = im2double(angle_temp,'indexed');
    angleResize = imresize(angle_temp,upsampling);
    
    center = round((sampling_size+1)/2);
    [psf_radius, ~] = size(angleResize);
    region = center+(-psf_radius/2:psf_radius/2-1);
    angle_1(region,region)=angleResize(:,:);
    angle_1 = ((angle_1./255)*2*pi);
else
    angle_1 = zeros(sampling_size, sampling_size);
    angle_temp = pmask;
    angleResize = imresize(angle_temp,upsampling);
    center = round((sampling_size+1)/2);
    [psf_radius, ~] = size(angleResize);
    region = center+(-psf_radius/2:psf_radius/2-1);
    angle_1(region,region)=angleResize(:,:);    
end
mask = angle_1;

A = 1/sqrt((pi/3)*(4+(rho_max^2-4)*sqrt(1-rho_max^2)));

gxx = sqrt(1./sqrt(1-rho.^2)).*(cos(phi).^2.*sqrt(1-rho.^2)+sin(phi).^2)*A; 
gyy = gxx.';
gxy = sqrt(1./sqrt(1-rho.^2)).*(sin(phi).*cos(phi).*(sqrt(1-rho.^2)-1))*A;
gyx = gxy.';
gxz = sqrt(1./sqrt(1-rho.^2)).*rho.*cos(phi)*A;
gyz = gxz.';

gxx(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gxy(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gxz(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gyx(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gyy(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;
gyz(xBFP.^2+yBFP.^2 > rho_max.^2) = 0;

gxPol = cat(3,gxx,gxy,gxz);
gyPol = cat(3,gyx,gyy,gyz);
[a,~,~] = size(gxPol);
bfp_radius = 165;
regionBFP = (round((a-1)/2)-(bfp_radius-1)/2):1:(round((a-1)/2)+(bfp_radius-1)/2);
basisgxPol = gxPol(regionBFP,regionBFP,:);
basisgyPol = gyPol(regionBFP,regionBFP,:);

switch BFPMod
    case 'xyPol'
        Gxx = fftshift(fft2(gxx));
        Gxy = fftshift(fft2(gxy));
        Gxz = fftshift(fft2(gxz));
        Gyx = fftshift(fft2(gyx));
        Gyy = fftshift(fft2(gyy));
        Gyz = fftshift(fft2(gyz));
       
    case 'pmask'
        Gxx = fftshift(fft2(gxx.*mask));
        Gxy = fftshift(fft2(gxy.*mask));
        Gxz = fftshift(fft2(gxz.*mask));
        Gyx = fftshift(fft2(gyx.*mask));
        Gyy = fftshift(fft2(gyy.*mask));
        Gyz = fftshift(fft2(gyz.*mask));

    case 'raPol'
        % Jones Vector
        J11 = cos(phi);
        J22 = cos(phi);
        J12 = sin(phi);
        J21 = -sin(phi);

        gxx_ra = gxx.*J11+gyx.*J12;
        gxy_ra = gxy.*J11+gyy.*J12;
        gxz_ra = gxz.*J11+gyz.*J12;
        gyx_ra = gxx.*J21+gyx.*J22;
        gyy_ra = gxy.*J21+gyy.*J22;
        gyz_ra = gxz.*J21+gyz.*J22;

        Gxx = fftshift(fft2(gxx_ra));
        Gxy = fftshift(fft2(gxy_ra));
        Gxz = fftshift(fft2(gxz_ra));
        Gyx = fftshift(fft2(gyx_ra));
        Gyy = fftshift(fft2(gyy_ra));
        Gyz = fftshift(fft2(gyz_ra)); 
end

GxPol = cat(3,Gxx,Gxy,Gxz);
GyPol = cat(3,Gyx,Gyy,Gyz);
[b,~,~] = size(GxPol);
region = (round((b-1)/2)-(image_size-1)/2):1:(round((b-1)/2)+(image_size-1)/2);
basisGxPol = GxPol(region,region,:);
basisGyPol = GyPol(region,region,:);

end

function f = circ(r)
f = (1 - sign(r-1))/2;
end
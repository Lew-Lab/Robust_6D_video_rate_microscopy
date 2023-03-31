function [bfpEx,bfpEy,imgEx,imgEy] = sim_Efield(Microscopy,pmask,BFPMod)

lambda = Microscopy.wavelength; n1 = Microscopy.n1; NA = Microscopy.NA;  n2 = Microscopy.n2;
z = Microscopy.z; z2 = Microscopy.z2; 
sampling_size = Microscopy.sampling_size;
upsampling = Microscopy.upsampling;
name = Microscopy.mask; rot=Microscopy.rot; xy_ind = Microscopy.xy_ind;

% generate the phase mask
if isnan(pmask)
    angle_temp = imread(name);
    angle_temp = rot90(angle_temp,rot);
    angle_1 = ones(sampling_size, sampling_size)*127;
    angle_temp = im2double(angle_temp,'indexed');
    angleResize = imresize(angle_temp,upsampling);
    
    center = round((sampling_size+1)/2);
    [psf_radius, ~] = size(angleResize);
    region = center+(-psf_radius/2:psf_radius/2-1);
    angle_1(region,region)=angleResize(:,:);
    angle_1 = ((angle_1/255)*2*pi)-pi;
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

[~,~,sizePmask] = size(mask);
if sizePmask==2
    pmaskx=mask(:,:,1);
    pmasky=mask(:,:,2);
elseif sizePmask==1
    pmaskx=mask;
    pmasky=mask;
end
if xy_ind==1
    pmaskx=rot90(pmask,3);
    pmasky=pmask;
end

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

% Exx: Ex contributed by mux
Exx = (1./sqrt(cos(theta1))).*(cos(phi).*Epx - sin(phi).*Esx).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
% Exy: Ex contributed by muy
Exy = (1./sqrt(cos(theta1))).*(cos(phi).*Epy - sin(phi).*Esy).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
% Exz: Ex contributed by muz
Exz = (1./sqrt(cos(theta1))).*(cos(phi).*Epz).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
% Eyx: Ey contributed by mux
Eyx = (1./sqrt(cos(theta1))).*(cos(phi).*Esx + sin(phi).*Epx).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
% Eyy: Ey contributed by muy
Eyy = (1./sqrt(cos(theta1))).*(cos(phi).*Esy + sin(phi).*Epy).*exp(1i*k1*z*cos(theta1)).*exp(1i*k2*z2*cos(theta2));
% Eyz: Ey contributed by muz
Eyz = (1./sqrt(cos(theta1))).*(sin(phi).*Epz).*exp(1i*k1*z*cos(theta1));

Exx(rho >= rho_max) = 0;
Exy(rho >= rho_max) = 0;
Exz(rho >= rho_max) = 0;
Eyx(rho >= rho_max) = 0;
Eyy(rho >= rho_max) = 0;
Eyz(rho >= rho_max) = 0;

bfpEx = cat(3,Exx,Exy,Exz);
bfpEy = cat(3,Eyx,Eyy,Eyz);

if ~ischar(BFPMod)
   error('MyComponent:incorrectType',...
       'Error. \nInput of BPF modulation must be a char, not a %s.',class(BFPMod))
end

switch BFPMod
    case 'xyPol'
        imgExx = fftshift(fft2(Exx));
        imgEyx = fftshift(fft2(Eyx));
        imgExy = fftshift(fft2(Exy));
        imgEyy = fftshift(fft2(Eyy));
        imgExz = fftshift(fft2(Exz));
        imgEyz = fftshift(fft2(Eyz));

        imgEx = cat(3,imgExx,imgExy,imgExz);
        imgEy = cat(3,imgEyx,imgEyy,imgEyz);

    case 'phase mask'
        imgExx = fftshift(fft2(Exx.*pmaskx));
        imgEyx = fftshift(fft2(Eyx.*pmasky));
        imgExy = fftshift(fft2(Exy.*pmaskx));
        imgEyy = fftshift(fft2(Eyy.*pmasky));
        imgExz = fftshift(fft2(Exz.*pmaskx));
        imgEyz = fftshift(fft2(Eyz.*pmasky));

        imgEx = cat(3,imgExx,imgExy,imgExz);
        imgEy = cat(3,imgEyx,imgEyy,imgEyz);

    case 'raPol'
        J11 = cos(phi);
        J22 = cos(phi);
        J12 = sin(phi);
        J21 = -sin(phi);

        bfpEx_out = zeros(size(bfpEx));
        bfpEy_out = zeros(size(bfpEy));
        imgEx = zeros(size(bfpEx));
        imgEy = zeros(size(bfpEy)); 
        for i = 1:3
            bfpEx_out(:,:,i) = bfpEx(:,:,i).*J11+bfpEy(:,:,i).*J12;
            bfpEy_out(:,:,i) = bfpEx(:,:,i).*J21+bfpEy(:,:,i).*J22;
            imgEx(:,:,i) = fftshift(fft2(bfpEx_out(:,:,i)));
            imgEy(:,:,i) = fftshift(fft2(bfpEy_out(:,:,i)));                
        end

    otherwise
        error('Please choose a valid BFP modulation')

end
% clear; clc;

% fname =  '2000_5_z1';
% 
% % offset = load([fname '_offset.mat']).offset;
% rawData = load([fname '_rawData.mat'],'-mat').rawData;
% % bg = load([fname '_bg.mat'],'-mat').b;
% locData = load([fname '_locData.mat'],'-mat').locData;
% 
% numFrame = size(rawData,3);
% 
% for f = 1:numFrame
%     img = squeeze(rawData(:,:,f));
%     if f == 1
%         imwrite(img,['rawData','.tif']);
%     else
%         imwrite(img,['rawData','.tif'],'writemode','append');
%     end
% end

z = 0.00000066662;
mu_x = 0.52568;
mu_y = -0.73108;
mu_z = 0.43495;
gamma = 1;
mxx = gamma*mu_x^2 + (1-gamma)/3;
myy = gamma*mu_y^2 + (1-gamma)/3;
mzz = gamma*mu_z^2 + (1-gamma)/3;
mxy = gamma*mu_x*mu_y;
mxz = gamma*mu_x*mu_z;
myz = gamma*mu_y*mu_z;

z_idx = round((z+1000e-9)/(6.6857e-08));
psf_layer = squeeze(dsf(:,:,z_idx,:,:));
img = zeros(8,101,101);
for ch = 1:8
    img(ch,:,:) = mxx*squeeze(psf_layer(ch,1,:,:)) + ...
        myy*squeeze(psf_layer(ch,2,:,:)) + ...
        mzz*squeeze(psf_layer(ch,3,:,:)) + ...
        mxy*squeeze(psf_layer(ch,4,:,:)) + ...
        mxz*squeeze(psf_layer(ch,5,:,:)) + ...
        myz*squeeze(psf_layer(ch,6,:,:));
end

img_cat = cat(2,squeeze(img(1,:,:)),squeeze(img(2,:,:)),squeeze(img(3,:,:)),...
    squeeze(img(4,:,:)),squeeze(img(5,:,:)),squeeze(img(6,:,:)),...
    squeeze(img(7,:,:)),squeeze(img(8,:,:)));

figure
imagesc(img_cat)
axis image
colormap bone
colorbar
title(['\mu_x=' num2str(mu_x) ', \mu_y=' num2str(mu_y) ', \mu_z=' num2str(mu_z) ', \gamma=' num2str(gamma)])
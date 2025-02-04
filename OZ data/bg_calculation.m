clear; clc;

fname =  '2000_5_z1';
offset = load([fname '_offset.mat']).offset;
rawData = load([fname '_rawData.mat'],'-mat').rawData;
bg = load([fname '_bg.mat'],'-mat').b;
numFrame = size(rawData,3);
data = zeros(numFrame,8,size(rawData,1),size(rawData,1));

for i = 1:numFrame
    dataTmp = double(rawData(:,:,i)) - offset;
    % camera brightness correction
    dataTmp(:,1:end/2) = dataTmp(:,1:end/2)*.49;
    dataTmp(:,end/2+1:end) = dataTmp(:,end/2+1:end)*.47;

    temp = zeros(8,size(dataTmp,1),size(dataTmp,1));
    for j = 1:8 
        temp(j,:,:) = dataTmp(:,size(dataTmp,1)*(j-1)+1:size(dataTmp,1)*j)';
        data(i,j,:,:) = temp(j,:,:);
    end
    data_cat(:,:,i) = cat(2,squeeze(temp(1,:,:)),squeeze(temp(2,:,:)),...
        squeeze(temp(3,:,:)),squeeze(temp(4,:,:)),...
        squeeze(temp(5,:,:)),squeeze(temp(6,:,:)),...
        squeeze(temp(7,:,:)),squeeze(temp(8,:,:)));

end

%%
data_avg = squeeze(sum(data,1))/size(data,1);
LPF = zeros(101,101);
LPF(50:52,50:52) = 1;
numFrame = size(data_cat,3);

figure(1)
subplot(3,1,2)
imagesc(data_cat(:,:,1))
colorbar
title('raw image')
axis image

for ch = 1:8
    BG = fftshift(fft2(squeeze(data_avg(ch,:,:))));
    BG = BG.*LPF;
    bg_lpf(ch,:,:) = abs(ifft2(ifftshift(BG)));
end

bg_cat = cat(2,squeeze(bg_lpf(1,:,:)),squeeze(bg_lpf(2,:,:)),...
    squeeze(bg_lpf(3,:,:)),squeeze(bg_lpf(4,:,:)),...
    squeeze(bg_lpf(5,:,:)),squeeze(bg_lpf(6,:,:)),...
    squeeze(bg_lpf(7,:,:)),squeeze(bg_lpf(8,:,:)));

bg = bg_lpf;
save('200_5_z1_bg_nonCat.mat',"bg");
bg = bg_cat;
save('200_5_z1_bg_cat.mat',"bg");

figure(1);
subplot(3,1,1)
imagesc(bg_cat)
title('bg')
axis image
colorbar

for f = 1:numFrame
    data_cat_bgfree(:,:,f) = squeeze(data_cat(:,:,f)) - bg_cat;
end

figure(1)
subplot(3,1,3)
imagesc
imagesc(data_cat_bgfree(:,:,1))
title('bgfree image')
colorbar
axis image
close all;

load(['data\',fName,'_offset.mat']);
load(['data\',fName,'.mat']);

img = 0;
k = 0;
for i = 50:100:size(rawData,3)
    img = img + double(rawData(:,:,i));
    k = k + 1;
end
img = img/k;

img = img - offset;
img(:,1:end/2) = img(:,1:end/2)*.49;
img(:,1+end/2:end) = img(:,1+end/2:end)*.47;
imgSz = size(img,1);

b = [];
for i = 1:8
    imgTmp = img(:,[1:imgSz]+imgSz*(i-1));
    bTmp = mean([imgTmp(1,:),imgTmp(end,:),imgTmp(:,1)',imgTmp(:,end)']);
    b = [b,bTmp];
end

figure();
subplot(3,1,1);
imagesc(img); axis image;
subplot(3,1,2);
imagesc(b); axis image;
caxis([min(img(:)),max(img(:))])
subplot(3,1,3);
imagesc(img-kron(b,ones(imgSz))); axis image;

save(['data\',fName,'_bg.mat'],'b');


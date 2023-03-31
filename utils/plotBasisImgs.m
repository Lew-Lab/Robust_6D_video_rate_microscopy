function plotBasisImgs(basisMatrix)

[a,~] = size(basisMatrix);
numPixel = sqrt(a/2);
% Bxx
temp=reshape(basisMatrix(:,1),numPixel,2*numPixel);
basisXXx=temp(:,1:numPixel);  basisXXy=temp(:,numPixel+1:2*numPixel);
basisXX=cat(1,basisXXx,basisXXy);
% Byy
temp=reshape(basisMatrix(:,2),numPixel,2*numPixel);
basisYYx=temp(:,1:numPixel);  basisYYy=temp(:,numPixel+1:2*numPixel);
basisYY=cat(1,basisYYx,basisYYy);
% Bzz
temp=reshape(basisMatrix(:,3),numPixel,2*numPixel);
basisZZx=temp(:,1:numPixel);  basisZZy=temp(:,numPixel+1:2*numPixel);
basisZZ=cat(1,basisZZx,basisZZy);
% Bxy
temp=reshape(basisMatrix(:,4),numPixel,2*numPixel);
basisXYx=temp(:,1:numPixel);  basisXYy=temp(:,numPixel+1:2*numPixel);
basisXY=cat(1,basisXYx,basisXYy);
% Bxz
temp=reshape(basisMatrix(:,5),numPixel,2*numPixel);
basisXZx=temp(:,1:numPixel);  basisXZy=temp(:,numPixel+1:2*numPixel);
basisXZ=cat(1,basisXZx,basisXZy);
% Byz
temp=reshape(basisMatrix(:,6),numPixel,2*numPixel);
basisYZx=temp(:,1:numPixel);  basisYZy=temp(:,numPixel+1:2*numPixel);
basisYZ=cat(1,basisYZx,basisYZy);

% Plot the Image Plane
figure('units','normalized','outerposition',[0 0 0.6 1]);
subplot(2,3,1);
imagesc(basisXX);
axis image; colorbar; axis off;
title("B_{xx}");

subplot(2,3,2);
imagesc(basisYY);
axis image; colorbar; axis off;
title("B_{yy}");

subplot(2,3,3);
imagesc(basisZZ);
axis image; colorbar; axis off;
title("B_{zz}");

subplot(2,3,4);
imagesc(basisXY);
axis image; colorbar; axis off;
title("B_{xy}");

subplot(2,3,5);
imagesc(basisXZ);
axis image; colorbar; axis off;
title("B_{xz}");

subplot(2,3,6);
imagesc(basisYZ);
axis image; colorbar; axis off;
title("B_{yz}");
end
function plotMVRBasisImgs(basisMatrix)

[a,~,~] = size(basisMatrix);
numPixel = a;

% Bxx
temp=basisMatrix(:,:,1);
basisXXr1 = temp(:,0*numPixel+1:numPixel);    basisXXr2 = temp(:,numPixel+1:2*numPixel);
basisXXr3 = temp(:,2*numPixel+1:3*numPixel);  basisXXr4 = temp(:,3*numPixel+1:4*numPixel);
basisXXa1 = temp(:,4*numPixel+1:5*numPixel);  basisXXa2 = temp(:,5*numPixel+1:6*numPixel);
basisXXa3 = temp(:,6*numPixel+1:7*numPixel);  basisXXa4 = temp(:,7*numPixel+1:8*numPixel);
% Byy
temp=basisMatrix(:,:,2);
basisYYr1 = temp(:,0*numPixel+1:numPixel);    basisYYr2 = temp(:,numPixel+1:2*numPixel);
basisYYr3 = temp(:,2*numPixel+1:3*numPixel);  basisYYr4 = temp(:,3*numPixel+1:4*numPixel);
basisYYa1 = temp(:,4*numPixel+1:5*numPixel);  basisYYa2 = temp(:,5*numPixel+1:6*numPixel);
basisYYa3 = temp(:,6*numPixel+1:7*numPixel);  basisYYa4 = temp(:,7*numPixel+1:8*numPixel);
% Bzz
temp=basisMatrix(:,:,3);
basisZZr1 = temp(:,0*numPixel+1:numPixel);    basisZZr2 = temp(:,numPixel+1:2*numPixel);
basisZZr3 = temp(:,2*numPixel+1:3*numPixel);  basisZZr4 = temp(:,3*numPixel+1:4*numPixel);
basisZZa1 = temp(:,4*numPixel+1:5*numPixel);  basisZZa2 = temp(:,5*numPixel+1:6*numPixel);
basisZZa3 = temp(:,6*numPixel+1:7*numPixel);  basisZZa4 = temp(:,7*numPixel+1:8*numPixel);
% Bxy
temp=basisMatrix(:,:,4);
basisXYr1 = temp(:,0*numPixel+1:numPixel);    basisXYr2 = temp(:,numPixel+1:2*numPixel);
basisXYr3 = temp(:,2*numPixel+1:3*numPixel);  basisXYr4 = temp(:,3*numPixel+1:4*numPixel);
basisXYa1 = temp(:,4*numPixel+1:5*numPixel);  basisXYa2 = temp(:,5*numPixel+1:6*numPixel);
basisXYa3 = temp(:,6*numPixel+1:7*numPixel);  basisXYa4 = temp(:,7*numPixel+1:8*numPixel);
% Bxz
temp=basisMatrix(:,:,5);
basisXZr1 = temp(:,0*numPixel+1:numPixel);    basisXZr2 = temp(:,numPixel+1:2*numPixel);
basisXZr3 = temp(:,2*numPixel+1:3*numPixel);  basisXZr4 = temp(:,3*numPixel+1:4*numPixel);
basisXZa1 = temp(:,4*numPixel+1:5*numPixel);  basisXZa2 = temp(:,5*numPixel+1:6*numPixel);
basisXZa3 = temp(:,6*numPixel+1:7*numPixel);  basisXZa4 = temp(:,7*numPixel+1:8*numPixel);
% Byz
temp=basisMatrix(:,:,6);
basisYZr1 = temp(:,0*numPixel+1:numPixel);    basisYZr2 = temp(:,numPixel+1:2*numPixel);
basisYZr3 = temp(:,2*numPixel+1:3*numPixel);  basisYZr4 = temp(:,3*numPixel+1:4*numPixel);
basisYZa1 = temp(:,4*numPixel+1:5*numPixel);  basisYZa2 = temp(:,5*numPixel+1:6*numPixel);
basisYZa3 = temp(:,6*numPixel+1:7*numPixel);  basisYZa4 = temp(:,7*numPixel+1:8*numPixel);

temp = hot(300); cMap(:,1) = min(temp(:,3),1); cMap(:,2) = min(temp(:,1)*0.50+temp(:,2)*0.50,1); cMap(:,3) = min(temp(:,1)*1,1);

figure('Position',[0 0 1600 1200])
ax1 = subplot(6,8,1); colormap(ax1,"hot");
imagesc(basisXXr1); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,2);
imagesc(basisXXr2); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,3);
imagesc(basisXXr3); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,4);
imagesc(basisXXr4); axis image; colorbar; axis off; colormap(ax1,"hot");
ax2 = subplot(6,8,5);
imagesc(basisXXa1); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,6);
imagesc(basisXXa2); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,7);
imagesc(basisXXa3); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,8);
imagesc(basisXXa4); axis image; colorbar; axis off; colormap(ax2,cMap);

ax1 = subplot(6,8,1*8+1);
imagesc(basisYYr1); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,1*8+2);
imagesc(basisYYr2); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,1*8+3);
imagesc(basisYYr3); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,1*8+4);
imagesc(basisYYr4); axis image; colorbar; axis off; colormap(ax1,"hot");
ax2 = subplot(6,8,1*8+5);
imagesc(basisYYa1); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,1*8+6);
imagesc(basisYYa2); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,1*8+7);
imagesc(basisYYa3); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,1*8+8);
imagesc(basisYYa4); axis image; colorbar; axis off; colormap(ax2,cMap);

ax1 = subplot(6,8,2*8+1);
imagesc(basisZZr1); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,2*8+2);
imagesc(basisZZr2); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,2*8+3);
imagesc(basisZZr3); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,2*8+4);
imagesc(basisZZr4); axis image; colorbar; axis off; colormap(ax1,"hot");
ax2 = subplot(6,8,2*8+5);
imagesc(basisZZa1); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,2*8+6);
imagesc(basisZZa2); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,2*8+7);
imagesc(basisZZa3); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,2*8+8);
imagesc(basisZZa4); axis image; colorbar; axis off; colormap(ax2,cMap);

ax1 = subplot(6,8,3*8+1);
imagesc(basisXYr1); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,3*8+2);
imagesc(basisXYr2); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,3*8+3);
imagesc(basisXYr3); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,3*8+4);
imagesc(basisXYr4); axis image; colorbar; axis off; colormap(ax1,"hot");
ax2 = subplot(6,8,3*8+5);
imagesc(basisXYa1); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,3*8+6);
imagesc(basisXYa2); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,3*8+7);
imagesc(basisXYa3); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,3*8+8);
imagesc(basisXYa4); axis image; colorbar; axis off; colormap(ax2,cMap);

ax1 = subplot(6,8,4*8+1);
imagesc(basisXZr1); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,4*8+2);
imagesc(basisXZr2); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,4*8+3);
imagesc(basisXZr3); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,4*8+4);
imagesc(basisXZr4); axis image; colorbar; axis off; colormap(ax1,"hot");
ax2 = subplot(6,8,4*8+5);
imagesc(basisXZa1); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,4*8+6);
imagesc(basisXZa2); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,4*8+7);
imagesc(basisXZa3); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,4*8+8);
imagesc(basisXZa4); axis image; colorbar; axis off; colormap(ax2,cMap);

ax1 = subplot(6,8,5*8+1);
imagesc(basisYZr1); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,5*8+2);
imagesc(basisYZr2); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,5*8+3);
imagesc(basisYZr3); axis image; colorbar; axis off; colormap(ax1,"hot");
ax1 = subplot(6,8,5*8+4);
imagesc(basisYZr4); axis image; colorbar; axis off; colormap(ax1,"hot");
ax2 = subplot(6,8,5*8+5);
imagesc(basisYZa1); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,5*8+6);
imagesc(basisYZa2); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,5*8+7);
imagesc(basisYZa3); axis image; colorbar; axis off; colormap(ax2,cMap);
ax2 = subplot(6,8,5*8+8);
imagesc(basisYZa4); axis image; colorbar; axis off; colormap(ax2,cMap);


h1 = annotation('textbox', [0.28 0.98 0 0], 'String', 'Radial Channel','FontSize',12, 'FitBoxToText', true,'EdgeColor','none');
h2 = annotation('textbox', [0.68 0.98 0 0], 'String', 'Azimuthal Channel','FontSize',12, 'FitBoxToText', true,'EdgeColor','none');
h3 = annotation('textbox', [0.1 0.89 0 0], 'String', 'B_{xx}','FontSize',12, 'FitBoxToText', true,'EdgeColor','none');
h4 = annotation('textbox', [0.1 0.75 0 0], 'String', 'B_{yy}','FontSize',12, 'FitBoxToText', true,'EdgeColor','none');
h3 = annotation('textbox', [0.1 0.61 0 0], 'String', 'B_{zz}','FontSize',12, 'FitBoxToText', true,'EdgeColor','none');
h4 = annotation('textbox', [0.1 0.47 0 0], 'String', 'B_{xy}','FontSize',12, 'FitBoxToText', true,'EdgeColor','none');
h3 = annotation('textbox', [0.1 0.33 0 0], 'String', 'B_{xz}','FontSize',12, 'FitBoxToText', true,'EdgeColor','none');
h4 = annotation('textbox', [0.1 0.19 0 0], 'String', 'B_{yz}','FontSize',12, 'FitBoxToText', true,'EdgeColor','none');
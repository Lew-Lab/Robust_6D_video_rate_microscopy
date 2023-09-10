clear; close all;

systemPar.zf = 1800e-9;
systemPar.zRoSE = [0,0,0];
systemPar.pxSize = 65*180/175*1e-9;
systemPar.lambda = 596e-9;
systemPar.PSFsz = 101;
systemPar.n1 = 1.515;
systemPar.n2 = 1.33;
systemPar.NA = 1.5;
systemPar.zList = (0:25:2000)*1e-9; 
systemPar.raRatio = 1;

Bstruct = MVRbasis_v4(systemPar);

%%
figure();
subplot(4,1,1);
imagesc(generateSphereImage(1,Bstruct,systemPar)); axis image;
title('\gamma = 1'); set(gca,'xtick',[],'ytick',[]);
subplot(4,1,2);
imagesc(generateSphereImage(.5,Bstruct,systemPar)); axis image;
title('\gamma = 0.5'); set(gca,'xtick',[],'ytick',[]);
subplot(4,1,3);
imagesc(generateSphereImage(.2,Bstruct,systemPar)); axis image;
title('\gamma = 0.2'); set(gca,'xtick',[],'ytick',[]);
subplot(4,1,4);
imagesc(generateSphereImage(0,Bstruct,systemPar)); axis image;
title('\gamma = 0'); set(gca,'xtick',[],'ytick',[]);

print('simSphere','-dpng');


%% compute image of 2-um sphere
function imgFull = generateSphereImage(gamma,Bstruct,systemPar)

pxSize = systemPar.pxSize*1e9;

imgFull = 0;
V = uniformSampleSphere(200);

for i = 1:size(V,1)
    x = V(i,1)*1e3;
    y = V(i,2)*1e3;
    z = V(i,3)*1e3+1e3;
    xIdx = round(x/pxSize);
    xRes = x - xIdx*pxSize;
    yIdx = round(y/pxSize);
    yRes = y - yIdx*pxSize;
    zIdx = round(z/25)+1;
    zRes = z - (zIdx-1)*25;

    mux = V(i,1);
    muy = V(i,2);
    muz = V(i,3);
    m = gamma*[mux^2,muy^2,muz^2,mux*muy,mux*muz,muz*muy]+(1-gamma)/3*[1,1,1,0,0,0];
    
    img = 0;
    for j = 1:6 
        % compute image including sub pixel shift
        % results would be similar even if we ignore subpixel shift since
        % sphere is much larger than one pixel
        img = img + m(j)*(Bstruct.Blist{zIdx}(:,:,j) + Bstruct.BgradxList{zIdx}(:,:,j)*xRes +...
             Bstruct.BgradyList{zIdx}(:,:,j)*yRes + Bstruct.BgradzList{zIdx}(:,:,j)*zRes);
        % img = img + m(j)*(Bstruct.Blist{zIdx}(:,:,j));
    end
    img = circshift(img,[yIdx,xIdx]); % I think you flipped x and y here,
    % the first element is shifting in the first dimension, which is
    % vertical direction of the matrix, and that's y...
    imgFull = imgFull + img;
end

end

%% unifromly sample the sphere
function V = uniformSampleSphere(N)

x = sin(acos(1-2*(0:(N-1))/(N-1))).*cos((0:(N-1))*2*pi*(1-2/(1+sqrt(5))));
y = sin(acos(1-2*(0:(N-1))/(N-1))).*sin((0:(N-1))*2*pi*(1-2/(1+sqrt(5))));
z = cos(acos(1-2*(0:(N-1))/(N-1)));
V = [x(:),y(:),z(:)];

end


function [B,xInd,yInd,zInd] = shiftBasis(x,y,z,pxSize,Blist,BgradxList,BgradyList,BgradzList,zStep,zList)

xInd = round(x/pxSize);
yInd = round(y/pxSize);
zInd = round((z-zList(1))/zStep)+1;

if abs(xInd) > size(Blist{1},1)/2
    xInd = ((size(Blist{1},1)/2)+.5)*sign(xInd);
end
if abs(yInd) > size(Blist{1},1)/2
    yInd = ((size(Blist{1},1)/2)+.5)*sign(yInd);
end
if zInd < 1
    zInd = 1;
end
if zInd > size(Blist,2)
    zInd = size(Blist,2);
end

psfSz = size(Blist{1},1);

Btmp(:,:,1:6) = Blist{zInd};
Btmp(:,:,7:9) = BgradxList{zInd}(:,:,1:3)*100;
Btmp(:,:,10:12) = BgradyList{zInd}(:,:,1:3)*100;
Btmp(:,:,13:15) = BgradzList{zInd}(:,:,1:3)*100;

for i = 1:8
    Btmp(:,psfSz*(i-1)+(1:psfSz),:) = circshift(Btmp(:,psfSz*(i-1)+(1:psfSz),:),[yInd,xInd]);
%     Btmp(:,psfSz*(i-1)+(1:psfSz),:) = circshift(Btmp(:,psfSz*(i-1)+(1:psfSz),:),yInd,1);
    if xInd > 0
        Btmp(:,psfSz*(i-1)+(1:xInd),:) = 0;
    elseif xInd <0
        Btmp(:,psfSz*(i-1)+(psfSz+xInd+1:psfSz),:) = 0;
    end
    if yInd > 0
        Btmp(1:yInd,psfSz*(i-1)+(1:psfSz),:) = 0;
    elseif yInd <0
        Btmp(psfSz+yInd+1:psfSz,psfSz*(i-1)+(1:psfSz),:) = 0;
    end
end

Bcrop = nan(psfSz-10,8*(psfSz-10),15);
for chNum = 1:8
    Bcrop(:,(chNum-1)*(psfSz-10)+(1:psfSz-10),:) = Btmp(6:psfSz-5,(chNum-1)*psfSz+(6:psfSz-5),:);
end   

B = reshape(Bcrop,[],15)';

end
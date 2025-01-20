function locData = driftCorrection(locData,c,r,Nz,Nxy,zThreshold,xyThreshold,ROIsz,vis)

locData = locData(:,locData(15,:)>500);
locData(2:4,:) = locData(2:4,:)*1e9 - c';
locData = locData(:,abs(locData(2,:))<ROIsz);
locData = locData(:,abs(locData(3,:))<ROIsz);

if strcmp(vis,'on')
    figure(); set(gcf,'position',[200,200,1500,800]);
    subplot(2,2,1);
    scatter(locData(2,:),locData(4,:),[],atand(locData(7,:)./locData(5,:)),'filled','sizedata',5);
    axis image;
end

%% correct z
for i = 1:ceil(locData(1,end)/Nz)
    frameIdx = logical((locData(1,:)<=i*Nz).*(locData(1,:)>(i-1)*Nz));
    locDataTmp = locData(:,frameIdx);
    locDataTmpZ = locDataTmp(:,locDataTmp(15,:)>1000);
    locDataTmpZ = locDataTmpZ(:,locDataTmpZ(2,:).^2+locDataTmpZ(3,:).^2>(r+50)^2);
    for iter = 1:10
        idxZ = abs(locDataTmpZ(4,:))<zThreshold;
        z0 = sum(locDataTmpZ(4,idxZ).*locDataTmpZ(15,idxZ))./sum(locDataTmpZ(15,idxZ));
        locDataTmpZ(4,:) = locDataTmpZ(4,:) - z0;
        locData(4,frameIdx) = locData(4,frameIdx) - z0;
    end
end

if strcmp(vis,'on')
    subplot(2,2,2);
    scatter(locData(2,:),locData(4,:),[],atand(locData(7,:)./locData(5,:)),'filled','sizedata',5);
    axis image;
end

%% correct xy 
for i = 1:ceil(locData(1,end)/Nxy)
    frameIdx = logical((locData(1,:)<=i*Nxy).*(locData(1,:)>(i-1)*Nxy));
    locDataTmp = locData(:,frameIdx);
    locDataTmpXY = locDataTmp(:,locDataTmp(15,:)>1000);
    locDataTmpXY = locDataTmpXY(:,abs(locDataTmpXY(4,:)-r)<xyThreshold);
    for iter = 1:10
        idxXY = logical( (locDataTmpXY(2,:).^2+locDataTmpXY(3,:).^2 > (r-50).^2) .* ...
            (locDataTmpXY(2,:).^2+locDataTmpXY(3,:).^2 < (r+50).^2) );
        c0 = fminsearch(@(cxy) sum((sqrt((locDataTmpXY(2,idxXY)-cxy(1)).^2+(locDataTmpXY(3,idxXY)-cxy(2)).^2)-r).^2.*...
            locDataTmpXY(15,idxXY)),[0,0]);        
        locDataTmpXY(2,:) = locDataTmpXY(2,:) - c0(1);
        locDataTmpXY(3,:) = locDataTmpXY(3,:) - c0(2);
        locData(2,frameIdx) = locData(2,frameIdx) - c0(1);
        locData(3,frameIdx) = locData(3,frameIdx) - c0(2);
    end
end

if strcmp(vis,'on')
    subplot(2,2,3);
    scatter(locData(2,:),locData(4,:),[],atand(locData(7,:)./locData(5,:)),'filled','sizedata',5);
    axis image; colormap hsv;
end

end
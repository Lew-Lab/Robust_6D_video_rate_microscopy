%clearvars -except fName fNameReg fNameSav;
clear;
fileFolder = 'C:\Users\cheng\OneDrive - Washington University in St. Louis\Desktop\Lew Lab\Robust-6D-video-rate-microscopy\experiment data\';
fileID = [5];
fNameRead = 'slb nfp 1000_crop';
fNameSave = 'ROI 3 NFP 1000';
fNameFull = [fileFolder,fNameRead,'.tif'];
locTmp = csvread(['09052023 reference','.csv'],1,0);
stacksNum = 9;

% find the ROI center of the 8 channels
locTmp = locTmp(:,2:3)/66.857;
for i = 1:8
    loc(i,:) = mean(locTmp(ceil(locTmp(:,1)/541)==i,:),1);
end
loc = round(loc);
loc(:,1) = loc(:,1) + 40;
loc(:,2) = loc(:,2) - 10;


m = 1;
imgSzHalf = 50;
frameNum = numel(imfinfo(fNameFull));
for n = 1:frameNum
    imgTmp = imread(fNameFull,n);
    img = [];
    for k = 1:8
        img = [img,imgTmp(loc(k,2)+[-imgSzHalf:imgSzHalf],loc(k,1)+[-imgSzHalf:imgSzHalf])];
    end
    rawData(:,:,m) = img;
    if m == 1
        imwrite(img,[fileFolder,fNameSave,'.tif']);
    else
        imwrite(img,[fileFolder,fNameSave,'.tif'],'writemode','append');
    end
    m = m + 1;
end

save([fileFolder,fNameSave,' raw.mat'],'rawData');

%%

offset = 0;
for i = 1:100
    imgTmp = imread([fileFolder,'offset.tif'],i);
    img = [];
    for j = 1:8
        img = [img,imgTmp(loc(j,2)+[-imgSzHalf:imgSzHalf],loc(j,1)+[-imgSzHalf:imgSzHalf])];
    end
    offset = offset + img;
end
offset = double(offset)/100;
xyRatio = 0.862275189291234;

temp = sum(rawData(:,:,1:100),3);
temp = temp - offset*100;

temp(:,1:end/2) = temp(:,1:end/2)*.49;
temp(:,end/2+1:end) = temp(:,end/2+1:end)*.47;

final = zeros(8,size(temp,1),size(temp,1));
for i = 1:8 
    final(i,:,:) = temp(:,101*(i-1)+1:101*i);
end

save([fNameSave '.mat'],'final');
imwrite(uint16(temp),[fileFolder,fNameSave,' final.tif']);
% imwrite(uint16(offset),[fileFolder,'slb 3.tif']);

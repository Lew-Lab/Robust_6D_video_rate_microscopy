clear

fileFolder = 'C:\Users\cheng\OneDrive - Washington University in St. Louis\Desktop\Lew Lab\Robust-6D-video-rate-microscopy\experiment data\';
fNameRead = 'fluorescent beads_crop';
fNameFull = [fileFolder,fNameRead,'.tif'];
fNameSave = 'registration_check';

% load the csv file and calculate the relative position
locTmp = csvread(['09052023 reference','.csv'],1,0);
locTmp = locTmp(:,2:3)/66.857;
for i = 1:8
    loc(i,:) = mean(locTmp(ceil(locTmp(:,1)/541)==i,:),1);
end
loc = round(loc);
loc_rela = loc - loc(1,:).*ones(size(loc));

chCenter_r = [1510 + loc_rela(1:5,1)'; 608 + loc_rela(1:5,2)'];

% load the tif and crop
m=1;
imgSzHalf = 70;
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


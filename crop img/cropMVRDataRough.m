clear; close all;
% tifpath = uigetdir();
tifpath = 'C:\Users\cheng\OneDrive - Washington University in St. Louis\Desktop\Lew Lab\Robust-6D-video-rate-microscopy\experiment data\';
ID_list = [8:10];
crop_size = 270;

chCenter_r = [1510,608;1449,1448;633,1357;676,587]';
chCenter_a = [1148,401;444,996;1085,1572;1733,1078]';

for ii = 1:length(ID_list)
%files = dir(fullfile([tifpath,num2str(ID_list(ii)),'\'],'*.tif'));
files = dir(fullfile([tifpath],['slb nfp 1000.tif']));

remFrameNew = 0;
for fInd = 1:numel(files)
    remFrame = remFrameNew;
    %fName = [[tifpath,num2str(ID_list(ii)),'\'],files(fInd).name];
    fName = [[tifpath],files(fInd).name];
    frameNum = numel(imfinfo(fName));
    if rem(frameNum,2)
        errordlg('stack contains odd number of frames')
        break;
    end
    for frameInd = 1:frameNum/2
        img1 = imread(fName,frameInd*2-1);
        img2 = imread(fName,frameInd*2);
        offset1 = sum(sum(img1(1:50,1:50)));
        offset2 = sum(sum(img2(1:50,1:50)));
        if offset1>offset2
            imgr = img1;
            imga = img2;
        else
            imgr = img2;
            imga = img1;
        end
        img = [];
        for i = 1:4
            img = [img,fliplr(imgr(round(chCenter_r(2,i))+[-crop_size:crop_size],round(chCenter_r(1,i))+[-crop_size:crop_size]))];
        end
        for i = 1:4
            img = [img,imga(round(chCenter_a(2,i))+[-crop_size:crop_size],round(chCenter_a(1,i))+[-crop_size:crop_size])];
        end
        if frameInd == 1
            imwrite(img,[fName(1:end-4),'_crop.tif']);
        else
            imwrite(img,[fName(1:end-4),'_crop.tif'],'writemode','append');
        end
    end
%     delete(fName);
end

end

%% generate mat form struct
load('B_ZScanned_zf1200_121.mat');
psfsize = systemPar.PSFsz; % set the size of the DSF
zsize = length(systemPar.zList);
cells = Bstruct.Blist;
dsf_raMVR = zeros(8,6,zsize,psfsize,psfsize);

for h = 1:zsize
    temp = cell2mat(cells(h));
    for o = 1:6
        for c = 1:8
            dsf_raMVR(c,o,h,:,:) = squeeze(temp(:,(c-1)*psfsize+1:c*psfsize,o));
        end
    end
end

dsf = dsf_raMVR;

save(['dsf_raMVR_psfSZ_' num2str(psfsize) '_zstack_' num2str(zsize)], 'dsf');

%% generate the video

videoname = 'dsf_MVR';
v = VideoWriter(videoname);
v.FrameRate = 1;
open(v);

for k = 1:zsize
    temp1 = cat(2,squeeze(dsf_raMVR(1,1,k,:,:)),squeeze(dsf_raMVR(2,1,k,:,:)),squeeze(dsf_raMVR(3,1,k,:,:)),squeeze(dsf_raMVR(4,1,k,:,:)),squeeze(dsf_raMVR(5,1,k,:,:)),squeeze(dsf_raMVR(6,1,k,:,:)),squeeze(dsf_raMVR(7,1,k,:,:)),squeeze(dsf_raMVR(8,1,k,:,:)));
    temp2 = cat(2,squeeze(dsf_raMVR(1,2,k,:,:)),squeeze(dsf_raMVR(2,2,k,:,:)),squeeze(dsf_raMVR(3,2,k,:,:)),squeeze(dsf_raMVR(4,2,k,:,:)),squeeze(dsf_raMVR(5,2,k,:,:)),squeeze(dsf_raMVR(6,2,k,:,:)),squeeze(dsf_raMVR(7,2,k,:,:)),squeeze(dsf_raMVR(8,2,k,:,:)));
    temp3 = cat(2,squeeze(dsf_raMVR(1,3,k,:,:)),squeeze(dsf_raMVR(2,3,k,:,:)),squeeze(dsf_raMVR(3,3,k,:,:)),squeeze(dsf_raMVR(4,3,k,:,:)),squeeze(dsf_raMVR(5,3,k,:,:)),squeeze(dsf_raMVR(6,3,k,:,:)),squeeze(dsf_raMVR(7,3,k,:,:)),squeeze(dsf_raMVR(8,3,k,:,:)));
    temp4 = cat(2,squeeze(dsf_raMVR(1,4,k,:,:)),squeeze(dsf_raMVR(2,4,k,:,:)),squeeze(dsf_raMVR(3,4,k,:,:)),squeeze(dsf_raMVR(4,4,k,:,:)),squeeze(dsf_raMVR(5,4,k,:,:)),squeeze(dsf_raMVR(6,4,k,:,:)),squeeze(dsf_raMVR(7,4,k,:,:)),squeeze(dsf_raMVR(8,4,k,:,:)));
    temp5 = cat(2,squeeze(dsf_raMVR(1,5,k,:,:)),squeeze(dsf_raMVR(2,5,k,:,:)),squeeze(dsf_raMVR(3,5,k,:,:)),squeeze(dsf_raMVR(4,5,k,:,:)),squeeze(dsf_raMVR(5,5,k,:,:)),squeeze(dsf_raMVR(6,5,k,:,:)),squeeze(dsf_raMVR(7,5,k,:,:)),squeeze(dsf_raMVR(8,5,k,:,:)));
    temp6 = cat(2,squeeze(dsf_raMVR(1,6,k,:,:)),squeeze(dsf_raMVR(2,6,k,:,:)),squeeze(dsf_raMVR(3,6,k,:,:)),squeeze(dsf_raMVR(4,6,k,:,:)),squeeze(dsf_raMVR(5,6,k,:,:)),squeeze(dsf_raMVR(6,6,k,:,:)),squeeze(dsf_raMVR(7,6,k,:,:)),squeeze(dsf_raMVR(8,6,k,:,:)));
    img = cat(1,temp1,temp2,temp3,temp4,temp5,temp6);
    imagesc(img);
    xlabel('second moments');
    ylabel('channel');
    axis image
    axis off
    F = getframe;
    writeVideo(v,F)
    
end

close(v);
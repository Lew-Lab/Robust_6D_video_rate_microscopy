%% generate mat form struct
psfsize = Microscopy.image_size; % set the size of the DSF
zsize = length(Microscopy.z2);
dsf_pixOL = zeros(2,6,zsize,psfsize,psfsize);

for h = 1:zsize
    Bxx = cell2mat(Bstruct.Bxx(h));
    Byy = cell2mat(Bstruct.Byy(h));
    Bzz = cell2mat(Bstruct.Bzz(h));
    Bxy = cell2mat(Bstruct.Bxy(h));
    Bxz = cell2mat(Bstruct.Bxz(h));
    Byz = cell2mat(Bstruct.Byz(h));
    for c = 1:2
        dsf_pixOL(c,1,h,:,:) = Bxx(:,45*(c-1)+1:45*(c-1)+45);
        dsf_pixOL(c,2,h,:,:) = Byy(:,45*(c-1)+1:45*(c-1)+45);
        dsf_pixOL(c,3,h,:,:) = Bzz(:,45*(c-1)+1:45*(c-1)+45);
        dsf_pixOL(c,4,h,:,:) = Bxy(:,45*(c-1)+1:45*(c-1)+45);
        dsf_pixOL(c,5,h,:,:) = Bxz(:,45*(c-1)+1:45*(c-1)+45);
        dsf_pixOL(c,6,h,:,:) = Byz(:,45*(c-1)+1:45*(c-1)+45);
    end
end

dsf = dsf_pixOL;
save(['dsf_pixOL_psfSZ_' num2str(psfsize) '_zstack_' num2str(zsize)], 'dsf');

%% generate the video

videoname = 'dsf_pixOL';
v = VideoWriter(videoname);
v.FrameRate = 1;
open(v);

for k = 1:zsize
    temp1 = cat(2,squeeze(dsf_pixOL(1,1,k,:,:)),squeeze(dsf_pixOL(2,1,k,:,:)));
    temp2 = cat(2,squeeze(dsf_pixOL(1,2,k,:,:)),squeeze(dsf_pixOL(2,2,k,:,:)));
    temp3 = cat(2,squeeze(dsf_pixOL(1,3,k,:,:)),squeeze(dsf_pixOL(2,3,k,:,:)));
    temp4 = cat(2,squeeze(dsf_pixOL(1,4,k,:,:)),squeeze(dsf_pixOL(2,4,k,:,:)));
    temp5 = cat(2,squeeze(dsf_pixOL(1,5,k,:,:)),squeeze(dsf_pixOL(2,5,k,:,:)));
    temp6 = cat(2,squeeze(dsf_pixOL(1,6,k,:,:)),squeeze(dsf_pixOL(2,6,k,:,:)));
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

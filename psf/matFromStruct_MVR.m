% generate mat form struct
file = 'B_zf1000_21PSFsz_2000zRange_200.571zStep_3px_sampling.mat';
load(file)
psfsize = systemPar.PSFsz; % Retrieve PSF size
zsize = length(systemPar.zList); % Retrieve PSF stack height

% Blist
cells = Bstruct.Blist;

%% Upsample the basis images if object grid space sampling != camera pixel size
scaling = floor(systemPar.pxSize/66.857e-09);
upsample_matrix = ones(scaling);
psfsize_new = psfsize * scaling;
    
% Returns array of DSFs with dimensions (Channels) x (Second Moments) x Z x X x Y
dsf = zeros(8,6,zsize,61,61); 

for h = 1:zsize
    temp = cell2mat(cells(h));
    for o = 1:6
        for c = 1:8
            dsf_slice = temp(:,(c-1)*psfsize+1:c*psfsize,o);
            dsf(c,o,h,:,:) = squeeze(dsf_slice);
        end
    end
end

% BsList
BsList = Bstruct.BsList;

% save dsf_raMVR as an MVR file with naming convention:
% MVR_zfxxx(the position of the NFP*1e9)_pixelsz/3pixelsz(sampling interval e.g. pixelsz is sampling at every pixel height, 3pixelsz is sampling at every 3 pixel heights)_zxxx(#z slices).mat
sampling = [num2str(floor(systemPar.pxSize/66.857e-09)) 'pixelsz'];
save(['MVR_zf' num2str(systemPar.zf*1e9) '_' sampling '_z' num2str(size(dsf, 3)) '.mat'],'dsf','Bstruct', 'systemPar');

%% generate the video

videoname = ['MVR_zf' num2str(systemPar.zf*1e9) '_' sampling '_z' num2str(size(dsf, 3))];
v = VideoWriter(videoname);
v.FrameRate = 1;
open(v);

for k = 1:zsize
     temp1 = cat(2,squeeze(dsf(1,1,k,:,:)),squeeze(dsf(2,1,k,:,:)),squeeze(dsf(3,1,k,:,:)),squeeze(dsf(4,1,k,:,:)),squeeze(dsf(5,1,k,:,:)),squeeze(dsf(6,1,k,:,:)),squeeze(dsf(7,1,k,:,:)),squeeze(dsf(8,1,k,:,:)));
     temp2 = cat(2,squeeze(dsf(1,2,k,:,:)),squeeze(dsf(2,2,k,:,:)),squeeze(dsf(3,2,k,:,:)),squeeze(dsf(4,2,k,:,:)),squeeze(dsf(5,2,k,:,:)),squeeze(dsf(6,2,k,:,:)),squeeze(dsf(7,2,k,:,:)),squeeze(dsf(8,2,k,:,:)));
     temp3 = cat(2,squeeze(dsf(1,3,k,:,:)),squeeze(dsf(2,3,k,:,:)),squeeze(dsf(3,3,k,:,:)),squeeze(dsf(4,3,k,:,:)),squeeze(dsf(5,3,k,:,:)),squeeze(dsf(6,3,k,:,:)),squeeze(dsf(7,3,k,:,:)),squeeze(dsf(8,3,k,:,:)));
     temp4 = cat(2,squeeze(dsf(1,4,k,:,:)),squeeze(dsf(2,4,k,:,:)),squeeze(dsf(3,4,k,:,:)),squeeze(dsf(4,4,k,:,:)),squeeze(dsf(5,4,k,:,:)),squeeze(dsf(6,4,k,:,:)),squeeze(dsf(7,4,k,:,:)),squeeze(dsf(8,4,k,:,:)));
     temp5 = cat(2,squeeze(dsf(1,5,k,:,:)),squeeze(dsf(2,5,k,:,:)),squeeze(dsf(3,5,k,:,:)),squeeze(dsf(4,5,k,:,:)),squeeze(dsf(5,5,k,:,:)),squeeze(dsf(6,5,k,:,:)),squeeze(dsf(7,5,k,:,:)),squeeze(dsf(8,5,k,:,:)));
     temp6 = cat(2,squeeze(dsf(1,6,k,:,:)),squeeze(dsf(2,6,k,:,:)),squeeze(dsf(3,6,k,:,:)),squeeze(dsf(4,6,k,:,:)),squeeze(dsf(5,6,k,:,:)),squeeze(dsf(6,6,k,:,:)),squeeze(dsf(7,6,k,:,:)),squeeze(dsf(8,6,k,:,:)));
    
     % Not sure why only the last five channels are being concatenated
     %{
    temp1 = cat(2,squeeze(dsf_raMVR(5,1,k,:,:)),squeeze(dsf_raMVR(6,1,k,:,:)),squeeze(dsf_raMVR(7,1,k,:,:)),squeeze(dsf_raMVR(8,1,k,:,:)));
    temp2 = cat(2,squeeze(dsf_raMVR(5,2,k,:,:)),squeeze(dsf_raMVR(6,2,k,:,:)),squeeze(dsf_raMVR(7,2,k,:,:)),squeeze(dsf_raMVR(8,2,k,:,:)));
    temp3 = cat(2,squeeze(dsf_raMVR(5,3,k,:,:)),squeeze(dsf_raMVR(6,3,k,:,:)),squeeze(dsf_raMVR(7,3,k,:,:)),squeeze(dsf_raMVR(8,3,k,:,:)));
    temp4 = cat(2,squeeze(dsf_raMVR(5,4,k,:,:)),squeeze(dsf_raMVR(6,4,k,:,:)),squeeze(dsf_raMVR(7,4,k,:,:)),squeeze(dsf_raMVR(8,4,k,:,:)));
    temp5 = cat(2,squeeze(dsf_raMVR(5,5,k,:,:)),squeeze(dsf_raMVR(6,5,k,:,:)),squeeze(dsf_raMVR(7,5,k,:,:)),squeeze(dsf_raMVR(8,5,k,:,:)));
    temp6 = cat(2,squeeze(dsf_raMVR(5,6,k,:,:)),squeeze(dsf_raMVR(6,6,k,:,:)),squeeze(dsf_raMVR(7,6,k,:,:)),squeeze(dsf_raMVR(8,6,k,:,:)));
    %}
    img = cat(1,temp1,temp2,temp3,temp4,temp5,temp6);
    imagesc(img);
    xlabel('second moments');
    ylabel('channel');
    axis image
    axis off
    colormap bone
    F = getframe();
    writeVideo(v,F)
    
end

close(v);

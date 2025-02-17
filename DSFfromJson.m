% MVR json config test
clc; close all; clear;
addpath matlabutils;

fname = "systemconfigs\systemParMVR.json"; % json file name
fid = fopen(fname); % open json file
raw = fread(fid,inf); % read json file 
str = char(raw');
fclose(fid); % once done reading, close file
systemPar = jsondecode(str); % decode json file contents into systemPar struct

% create zList from systemPar zRange and zStep
systemPar.zList = round(systemPar.zRange(1):systemPar.zStep:systemPar.zRange(2),11); 

%% Generate: Basis images and basis image gradients 

wb = waitbar(0,'basis generation step (1/3) 0%');
Bstruct = MVRbasis_v4(systemPar,wb);
close(wb);

%% Generate DSFs from basis images

cells = Bstruct.Blist;

zsize = length(systemPar.zList); % Retrieve PSF stack height

psfsize = systemPar.PSFsz;

% Returns array of DSFs with dimensions (Channels) x (Second Moments) x Z x
% Y x X (?, ask Jasmine about this)
dsf = zeros(8,6,zsize,psfsize,psfsize); 

for h = 1:zsize
    % basis images for particular z height; 
    % dimensions (psfsize) x (psfsize * 8) x 6
    temp = cell2mat(cells(h)); 

    for c = 1:8
        for o = 1:6
            dsf_slice = temp(:,(c-1)*psfsize+1:c*psfsize,o);
            dsf_temp = squeeze(dsf_slice);
            dsf(c,o,h,:,:) = dsf_temp;
        end

    end
end


% save dsf_raMVR as an MVR file with naming convention:
% MVR_zfxxx(the position of the NFP*1e9)_pixelsz/3pixelsz(sampling interval e.g. pixelsz is sampling at every pixel height, 3pixelsz is sampling at every 3 pixel heights)_zxxx(#z slices).mat
sampling = [num2str(floor(systemPar.pxSize/66.857e-09)) 'pixelsz'];
savename = ['psf\MVR_zf' num2str(systemPar.zf*1e9) '_' num2str(systemPar.PSFsz) 'PSFsz_' sampling '_z' num2str(size(dsf, 3)) '.mat'];
save(savename,'dsf','Bstruct', 'systemPar');

%% generate the video

videoname = savename;
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

    % retrieve DSF slices at top, middle, and bottom of the generated
    % stack.
    if k == 1
        dsf_bottom = img;
    elseif k == floor(zsize/2)
        dsf_mid = img;
    elseif k == zsize
        dsf_top = img;
    end
      
    xlabel('second moments');
    ylabel('channel');
    axis image
    axis off
    colormap bone
    F = getframe();
    writeVideo(v,F)
    
end

close(v);


%% Save basis images & gradients + system parameters

savename = ['B_zf' num2str(systemPar.zf*1e9) '_' num2str(systemPar.PSFsz) ...
    'PSFsz_' num2str(systemPar.zRange(2)*1e9) 'zRange_' num2str(systemPar.zStep*1e9) 'zStep_' num2str(floor(systemPar.pxSize/66.857e-09)) 'px_sampling'];
save(['psf\' savename '.mat'],'Bstruct','systemPar');
 
%% DSF intensity profile analysis (bottom, top, middle)

intensity_profile_dir = 'psf\psf intensity profiles\';
mkdir([intensity_profile_dir savename])
% Generate intensity profile for the DSFs in each channel and second moment

for c = 1:8
    for s = 1:6
        cs_bottom = dsf_bottom((s-1)*psfsize+1:s*psfsize,(c-1)*psfsize+1:c*psfsize);
        cs_mid = dsf_mid((s-1)*psfsize+1:s*psfsize,(c-1)*psfsize+1:c*psfsize);
        cs_top = dsf_top((s-1)*psfsize+1:s*psfsize,(c-1)*psfsize+1:c*psfsize);
        
        cs_cell_temp = cell(1, 3); 
        cs_cell_temp{1} = cs_bottom; cs_cell_temp{2} = cs_mid; cs_cell_temp{3} = cs_top;

        for k = 1:size(cs_cell_temp, 2)
            % find global max in channel
            cs_temp = cs_cell_temp{k};
            [max_val, max_idx_temp] = max(abs(cs_temp),[],"all");
            [row, col] = ind2sub([psfsize psfsize], max_idx_temp);
            
            secname = ['xx';'yy';'zz';'xy';'xz';'yz'];
            figname = ['channel ' num2str(c) ' m_{' num2str(secname(s,:)) '}'];
            zpos = {'z_{bottom}'; 'z_{mid}'; 'z_{top}'};

            f = figure(2); set(gcf, 'Position', [0 0 3000 700]);
            subplot(1,3,1); hold on;
            imagesc(cs_temp);
            line([1 psfsize], [row row], 'Color', 'b', 'LineWidth', 2);
            line([col col], [1 psfsize],'Color', 'r', 'LineWidth', 2);
            axis image ij; colormap bone; colorbar;
            title([figname ', ' zpos{k}]);

            subplot(1,3,2)
            improfile(cs_temp, [1 psfsize], [row row]);
            xlim([1 psfsize]); xlabel("Distance along profile (px)")
            axis square;
            title([figname ', horiz, ' zpos{k}])

            subplot(1,3,3)
            [a,b,colr] = improfile(cs_temp, [col col], [1 psfsize]);
            plot(b,colr,'Color', 'r'); xlabel('Distance along profile (px)')
            xlim([1 psfsize])
            axis square;
            title([figname ', vert, ' zpos{k}])

            savefig(f, [intensity_profile_dir savename '\' figname ', ' zpos{k} ' figure.fig']);
            saveas(f, [intensity_profile_dir savename '\' figname  ', ' zpos{k} ' img.png']);
        end
    end
end
            




            

        
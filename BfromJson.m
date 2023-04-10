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

%% Save basis images & gradients + system parameters

save(['psf' filesep 'B_ZScanned_zf' num2str(systemPar.zf*1e9) '_PSFsz' num2str(systemPar.PSFsz) '_zList' num2str(length(systemPar.zList))],'Bstruct','systemPar','-v7.3');


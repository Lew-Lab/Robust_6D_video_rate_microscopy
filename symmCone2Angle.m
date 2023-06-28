%%
addpath utils;
obj = load(['objects/sphere 2.mat']);
object = obj.object;

%% Load and format second moment information for input

% Format of estimated data: 6 x Z x X x Y 
% (6 rows are m_xx, m_yy, m_zz, m_xy, m_xz, m_yz)
% each Z x X x Y would then be m_xx, m_yy, m_zz, m_xy, m_xz, m_yz at
% each voxel
m_xx = squeeze(object(1,:,:,:));
m_yy = squeeze(object(2,:,:,:));
m_zz = squeeze(object(3,:,:,:));

s = m_xx+m_yy+m_zz; s = reshape(s, [], 1);

m_xy = squeeze(object(4,:,:,:));
m_xz = squeeze(object(5,:,:,:));
m_yz = squeeze(object(6,:,:,:));

m2 = [reshape(m_xx, [], 1), reshape(m_yy, [], 1), reshape(m_zz, [], 1), reshape(m_xy, [], 1), reshape(m_xz, [], 1), reshape(m_yz, [], 1)];

%% retrieve indices of nonzero second moments
[z,x,y] = ind2sub(size(m_xx), find(s > 1e-2)); % find indices of localizations that exist
m2 = m2((s >1e-2),:); % remove second moments of localizations that don't exist i.e. localizations with all zero second moments
%% Load imaging system basis images

b = load(['psf\B_ZScanned_zf1000_PSFsz_61.mat']);
Bstruct = b.Bstruct;

BsList = Bstruct.BsList;
B_HaList = Bstruct.B_HaList;
sumNormList = Bstruct.sumNormList;

% B is known from the imaging system and parameters
% background is also known?
m1 = zeros(length(m2),4);
for loc = 1:length(m2)
    m2Tmp = m2(loc, 1:6); % second moments associated with localization
    signal = sum(m2Tmp(1:3)); % where mTmp = s*(all second moments). Sum the first three (s*mii) to recover s

    secM = m2Tmp(1:6); % normalized second moments associated with localization

    [muxTmp, muyTmp, muzTmp, gammaTmp, x0] = secondM2SymmConeWeighted(BsList{z(loc)},B_HaList{z(loc)},...
                    sumNormList{z(loc)},secM/signal,signal,1e-12);
    %[muxTmp, muyTmp, muzTmp, gammaTmp, x0] = secondM2SymmConeWeighted(BsList, B_HaList, sumNormList, secM, signal, zeros(1,6));
    m1(loc,1:4) = [muxTmp muyTmp abs(muzTmp) gammaTmp]; % fix mu_z to be positive so calculated theta is positive
end

%% recovering angles and wobble:
% mu_x = sin(theta)cos(phi)
% mu_y = sin(theta)sin(phi)
% mu_z = cos(theta) -> theta = arccos(mu_z)
mux = m1(:,1);
muy = m1(:,2);
muz = m1(:,3);
gamma = m1(:,4);

% theta and phi in degrees
theta = acosd(muz);
phi = atand(muy./mux);

% convert rotational mobility (gamma) to Omega
% gamma = 1-(3*Omega)/4*pi + Omega^2/(8*pi^2)
% quadratic equation: Omega^2/(8*pi^2) - (3*Omega)/4*pi + (1-gamma)

omega = zeros(length(gamma),2);
for loc = 1:length(gamma)
    omega(loc,:) = roots([1/(8*pi^2) -3/(4*pi) (1-gamma(loc))]);
end

%% plot first moment data: 3D

figure(1);
scatter3(x, y, z,[], phi,'filled','markerfacealpha',.4)
xlabel('x (voxel length)'); ylabel('y (voxel length)'); zlabel('z (voxel length)');
title('\phi: 3D diagonal view');
colormap hsv; colorbar;

figure(2);
scatter3(x, y, z,[], theta,'filled','markerfacealpha',.4)
xlabel('x (voxel length)'); ylabel('y (voxel length)'); zlabel('z (voxel length)');
title('\theta: 3D diagonal view');
colormap cool; colorbar;

%% plot first moment data: cross sections
figure(3);
scatter(x,y,[], phi,'filled','markerfacealpha',.4);
xlabel('x (voxel length)'); ylabel('y (voxel length)'); axis square;
title("\phi"); colormap hsv; colorbar;

figure(4);
scatter(x,y,[],theta,'filled','markerfacealpha',.4);
xlabel('x (voxel length)'); ylabel('y (voxel length)'); axis square;
title("\theta"); colormap cool; colorbar;


figure(5);
scatter(x,z,[], phi,'filled','markerfacealpha',.4);
xlabel('x (voxel length)'); ylabel('z (voxel length)'); axis square;
title("\phi"); colormap hsv; colorbar;


figure(6);
scatter(x,z,[],theta,'filled','markerfacealpha',.4);
xlabel('x (voxel length)'); ylabel('z (voxel length)'); axis square;
title("\theta"); colormap cool; colorbar;


figure(7);
scatter(y,z,[], phi,'filled','markerfacealpha',.4);
xlabel('y (voxel length)'); ylabel('z (voxel length)'); axis square;
title("\phi"); colormap hsv; colorbar;


figure(8);
scatter(y,z,[],theta,'filled','markerfacealpha',.4);
xlabel('y (voxel length)'); ylabel('z (voxel length)'); axis square;
title("\theta"); colormap cool; colorbar;

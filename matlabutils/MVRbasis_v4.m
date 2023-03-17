function Bstruct = MVRbasis_v4(systemPar,varargin)

% varargin - waitbar

pxSize = systemPar.pxSize; % in meter
lambda = systemPar.lambda; % in meter
PSFsz = systemPar.PSFsz; % # of pixels
n1 = systemPar.n1; % imaging medium ri
n2 = systemPar.n2; % sample ri
NA = systemPar.NA;
zf = systemPar.zf; % focal plane position (positive above the sample)
zList = systemPar.zList; % basis dictionary z position list
zRoSE = systemPar.zRoSE; % 1x3 vector [min, approx focal, max] to calculate the average PSF for RoSE
raRatio = systemPar.raRatio; % radial channel * raRatio

if ~isempty(varargin)
    wb = waitbar(0,varargin{1},'basis generation step (1/3) 0%');
end

for i = 1:length(zList)
    z = zList(i);
    [~,B,Bgradx,Bgrady,Bgradz,Bs,B_Ha,sumNorm] = computeBasis(n1, n2, NA, z, zf, pxSize, lambda, PSFsz, raRatio);
    Blist{i} = B;
    BgradxList{i} = Bgradx;
    BgradyList{i} = -Bgrady;
    BgradzList{i} = Bgradz;
    BsList{i} = Bs;
    B_HaList{i} = B_Ha;
    sumNormList{i} = sumNorm;
    if exist('wb','var')
        wb = waitbar(i/length(zList),wb,['basis generation step (1/3) ',num2str(round(i/length(zList)*100)),'%']);
    end
end

% prepare FPSF structure for RoSE
slopexz = 1.37 - NA/5; % approximate, probably need futher correction
%% azimuthal channel

zStep_a = pxSize/slopexz;
zList_a_tmp1 = zRoSE(2): -zStep_a : zRoSE(1);
zList_a_tmp2 = zRoSE(2): zStep_a : zRoSE(3);
zList_a = [fliplr(zList_a_tmp1),zList_a_tmp2(2:end)];
centerIdx = numel(zList_a_tmp1);
B_a = zeros(PSFsz);
Bgradx_a = zeros(PSFsz);
Bgrady_a = zeros(PSFsz);
for i = 1:length(zList_a)
    z = zList_a(i);
    [~,B,Bgradx,Bgrady] = computeBasis(n1, n2, NA, z, zf, pxSize, lambda, PSFsz, raRatio);
    
    % create azimuthally polarized basis fields from B
    B_a = B_a + circshift(sum(B(:,PSFsz*4+(1:PSFsz),1:3),3),i-centerIdx);
    Bgradx_a = Bgradx_a + circshift(sum(Bgradx(:,PSFsz*4+(1:PSFsz),1:3),3),i-centerIdx);
    Bgrady_a = Bgrady_a + circshift(sum(Bgrady(:,PSFsz*4+(1:PSFsz),1:3),3),i-centerIdx);
    if exist('wb','var')
        wb = waitbar(i/length(zList_a),wb,['basis generation step (2/3) ',num2str(round(i/length(zList_a)*100)),'%']);
    end
end
Bscaling = sum(B_a(:)); % scaling / normalization factor
B_a = B_a/Bscaling;
Bgradx_a = Bgradx_a/Bscaling;
Bgrady_a = Bgrady_a/Bscaling;

FPSF_a.FXX = (fft2((fftshift(B_a))));
FPSF_a.FXXdx = (fft2((fftshift(10^2*Bgradx_a))));
FPSF_a.FXXdy = -(fft2((fftshift(10^2*Bgrady_a)))); % negative to match RoSE-O coordinates

%% radial channel

zStep_r = pxSize/slopexz*sqrt(2);
zList_r_tmp1 = zRoSE(2): -zStep_r : zRoSE(1);
zList_r_tmp2 = zRoSE(2): zStep_r : zRoSE(3);
zList_r = [fliplr(zList_r_tmp1),zList_r_tmp2(2:end)];
centerIdx = numel(zList_r_tmp1);
B_r = zeros(PSFsz);
Bgradx_r = zeros(PSFsz);
Bgrady_r = zeros(PSFsz);
for i = 1:length(zList_r)
    z = zList_r(i);
    [~,B,Bgradx,Bgrady] = computeBasis(n1, n2, NA, z, zf, pxSize, lambda, PSFsz, raRatio);
    
    % create radially polarized basis fields from B
    B_r = B_r + circshift(sum(B(:,(1:PSFsz),1:3),3),[i-centerIdx,i-centerIdx]);
    Bgradx_r = Bgradx_r + circshift(sum(Bgradx(:,(1:PSFsz),1:3),3),[i-centerIdx,i-centerIdx]);
    Bgrady_r = Bgrady_r + circshift(sum(Bgrady(:,(1:PSFsz),1:3),3),[i-centerIdx,i-centerIdx]);
    if exist('wb','var')
        wb = waitbar(i/length(zList_r),wb,['basis generation step (3/3) ',num2str(round(i/length(zList_r)*100)),'%']);
    end
end
Bscaling = sum(B_r(:));
B_r = B_r/Bscaling;
Bgradx_r = Bgradx_r/Bscaling;
Bgrady_r = Bgrady_r/Bscaling;

FPSF_r.FXX = (fft2((fftshift(B_r))));
FPSF_r.FXXdx = (fft2((fftshift(10^2*Bgradx_r))));
FPSF_r.FXXdy = -(fft2((fftshift(10^2*Bgrady_r)))); % negative to match RoSE-O coordinates

% save struct
Bstruct.Blist = Blist;
Bstruct.BgradxList = BgradxList;
Bstruct.BgradyList = BgradyList;
Bstruct.BgradzList = BgradzList;
Bstruct.FPSF_r = FPSF_r;
Bstruct.FPSF_a = FPSF_a;
Bstruct.BsList = BsList;
Bstruct.B_HaList = B_HaList;
Bstruct.sumNormList = sumNormList;

end

%%
function [E,B,Bgradx,Bgrady,Bgradz,Bs,B_Ha,sumNorm] = computeBasis(n1, n2, NA, z, zf, pxSize, lambda, PSFsz, raRatio)

% zf - position of nfp above the ri interface
% n1 - immersion oil r.i.
% n2 - sample r.i.

rm = NA/n1;
bfpRadius = 40; 
bfpN = ceil(bfpRadius*2/rm);
if rem(bfpN,2) == 0
    bfpN = bfpN + 1;
end
[x,y] = meshgrid(linspace(-1,1,bfpN),linspace(-1,1,bfpN));
[t,r] = cart2pol(x,y);

N = ceil(lambda/pxSize/NA*bfpRadius);
if rem(N,2) == 0
    N = N + 1;
end

t_zeroPad = zeros(N);
t_zeroPad((N-bfpN)/2+(1:bfpN),(N-bfpN)/2+(1:bfpN)) = t;

J11 = cos(t_zeroPad);
J22 = cos(t_zeroPad);
J12 = sin(t_zeroPad);
J21 = -sin(t_zeroPad);

xChPart = zeros(N);
xChPart((N-bfpN)/2+(1:bfpN),(N-bfpN)/2+(1:bfpN)) = (x>0).*(y>0);
yChPart = zeros(N);
yChPart((N-bfpN)/2+(1:bfpN),(N-bfpN)/2+(1:bfpN)) = (t>pi/4).*(t<pi*3/4);

[X,~] = meshgrid(1e-9*(1:N),1e-9*(1:N));
maskdx = exp(1j*2*pi*X/(N*pxSize));
maskdy = rot90(maskdx);
maskdzTmp = exp(1j*2*pi*n2/lambda*1e-9*cos(asin(n1/n2*r)));
maskdz = zeros(N);
maskdz((N-bfpN)/2+(1:bfpN),(N-bfpN)/2+(1:bfpN)) = maskdzTmp;

%%
dx = n1*pxSize;
[eta,xi] = meshgrid(linspace(-1/(2*dx),1/(2*dx),N),linspace(-1/(2*dx),1/(2*dx),N)); % image plane coordinate space?

xBFP = lambda*eta;
yBFP = lambda*xi;
[phi,rho] = cart2pol(xBFP,yBFP);
rho_max = NA/n1;

k1 = n1*(2*pi/lambda);
k2 = n2*(2*pi/lambda);
theta1 = asin(rho);
theta2 = asin((n1/n2)*sin(theta1));

a = cos(theta2)./cos(theta1);

%Fresnel coefficients
ts = 2*n1./(n1+n2*a);
tpxy = 2*n1./(n1+n2./a);
tpz = 2*n1./(n2+n1*a);

Esx = -sin(phi).*ts;
Esy = cos(phi).*ts;
Epx = cos(phi).*cos(theta1).*tpxy;
Epy = sin(phi).*cos(theta1).*tpxy;
Epz = -sin(theta1).*(n1/n2).*tpz;

% Exx - Ex contributed by mux
Exx = (1./sqrt(cos(theta1))).*(cos(phi).*Epx - sin(phi).*Esx).*exp(1i*k2*z*cos(theta2)).*exp(-1i*k1*zf*cos(theta1));
Exy = (1./sqrt(cos(theta1))).*(cos(phi).*Epy - sin(phi).*Esy).*exp(1i*k2*z*cos(theta2)).*exp(-1i*k1*zf*cos(theta1));
Exz = (1./sqrt(cos(theta1))).*(cos(phi).*Epz).*exp(1i*k2*z*cos(theta2)).*exp(-1i*k1*zf*cos(theta1));
Eyx = (1./sqrt(cos(theta1))).*(cos(phi).*Esx + sin(phi).*Epx).*exp(1i*k2*z*cos(theta2)).*exp(-1i*k1*zf*cos(theta1));
Eyy = (1./sqrt(cos(theta1))).*(cos(phi).*Esy + sin(phi).*Epy).*exp(1i*k2*z*cos(theta2)).*exp(-1i*k1*zf*cos(theta1));
Eyz = (1./sqrt(cos(theta1))).*(sin(phi).*Epz).*exp(1i*k2*z*cos(theta2)).*exp(-1i*k1*zf*cos(theta1));

Exx(rho >= rho_max) = 0;
Exy(rho >= rho_max) = 0;
Exz(rho >= rho_max) = 0;
Eyx(rho >= rho_max) = 0;
Eyy(rho >= rho_max) = 0;
Eyz(rho >= rho_max) = 0;

Ex_bfp = cat(3,Exx,Exy,Exz); % concatenate along 3D dimension
Ey_bfp = cat(3,Eyx,Eyy,Eyz); % concatenate along 3D dimension

%%
PSFregion = N/2+[-PSFsz/2+1:PSFsz/2];

for j = 1:4
    for i = 1:3
        if j == 1
            Ex_bfpTmp(:,:,i) = Ex_bfp(:,:,i);
            Ey_bfpTmp(:,:,i) = Ey_bfp(:,:,i);
        elseif j == 2 % gradient w.r.t x
            Ex_bfpTmp(:,:,i) = Ex_bfp(:,:,i).*maskdx;
            Ey_bfpTmp(:,:,i) = Ey_bfp(:,:,i).*maskdx;
        elseif j == 3 % gradient w.r.t y
            Ex_bfpTmp(:,:,i) = Ex_bfp(:,:,i).*maskdy;
            Ey_bfpTmp(:,:,i) = Ey_bfp(:,:,i).*maskdy;
        else % gradient w.r.t z
            Ex_bfpTmp(:,:,i) = Ex_bfp(:,:,i).*maskdz;
            Ey_bfpTmp(:,:,i) = Ey_bfp(:,:,i).*maskdz;
        end
        Ex_bfp_J(:,:,i) = Ex_bfpTmp(:,:,i).*J11+Ey_bfpTmp(:,:,i).*J12;
        Ey_bfp_J(:,:,i) = Ex_bfpTmp(:,:,i).*J21+Ey_bfpTmp(:,:,i).*J22;
        Ex(:,:,i) = [fftshift(fft2(Ex_bfp_J(:,:,i).*xChPart)),fftshift(fft2(Ex_bfp_J(:,:,i).*rot90(xChPart))),...
            fftshift(fft2(Ex_bfp_J(:,:,i).*rot90(xChPart,2))),fftshift(fft2(Ex_bfp_J(:,:,i).*rot90(xChPart,3)))]*sqrt(raRatio);
        Ey(:,:,i) = [fftshift(fft2(Ey_bfp_J(:,:,i).*yChPart)),fftshift(fft2(Ey_bfp_J(:,:,i).*rot90(yChPart))),...
            fftshift(fft2(Ey_bfp_J(:,:,i).*rot90(yChPart,2))),fftshift(fft2(Ey_bfp_J(:,:,i).*rot90(yChPart,3)))];
    end
    Elist{j} = [Ex(PSFregion,[PSFregion,PSFregion+N,PSFregion+N*2,PSFregion+N*3],:),...
        Ey(PSFregion,[PSFregion,PSFregion+N,PSFregion+N*2,PSFregion+N*3],:)];
end
E = Elist{1};

for j = 1:4 % compute B (1) and gradients of B w.r.t x (2), y (3), and z (4)
    E_tmp = Elist{j};
    for i = 1:3 % compute Bxx, Byy, Bzz
        B_tmp(:,:,i) = abs(E_tmp(:,:,i)).^2; 
    end
    % compute Bxy, Bxz, Byz
    B_tmp(:,:,4) = 2*real(E_tmp(:,:,1).*conj(E_tmp(:,:,2)));
    B_tmp(:,:,5) = 2*real(E_tmp(:,:,1).*conj(E_tmp(:,:,3)));
    B_tmp(:,:,6) = 2*real(E_tmp(:,:,2).*conj(E_tmp(:,:,3)));
    B_tmp = B_tmp/sum(sum(B_tmp(:,:,1))); % normalization factor?
    Blist{j} = B_tmp;
end
B = Blist{1};

% calculate B gradients
Bgradx = Blist{2}-Blist{1};
Bgrady = Blist{3}-Blist{1};
Bgradz = Blist{4}-Blist{1};

if PSFsz > 31
    B0 = [];
    for i = 0:7
        B0 = [B0,B((PSFsz+1)/2+[-15:15],PSFsz*i+(PSFsz+1)/2+[-15:15],:)];
    end
else
    B0 = B;
end
% prepare first moment projection
Bs.XX = B0(:,:,1);
Bs.YY = B0(:,:,2);
Bs.ZZ = B0(:,:,3);
Bs.XY = B0(:,:,4);
Bs.XZ = B0(:,:,5);
Bs.YZ = B0(:,:,6);

sumNorm = sum(sum(Bs.XX));

% Hadamard (aka element-wise) products
B_Ha.aa = (Bs.XX) .* (Bs.XX);
B_Ha.ab = (Bs.XX) .* (Bs.YY);
B_Ha.ac = (Bs.XX) .* (Bs.ZZ);
B_Ha.ad = (Bs.XX) .* (Bs.XY);
B_Ha.ae = (Bs.XX) .* (Bs.XZ);
B_Ha.af = (Bs.XX) .* (Bs.YZ);

B_Ha.bb = (Bs.YY) .* (Bs.YY);
B_Ha.bc = (Bs.YY) .* (Bs.ZZ);
B_Ha.bd = (Bs.YY) .* (Bs.XY);
B_Ha.be = (Bs.YY) .* (Bs.XZ);
B_Ha.bf = (Bs.YY) .* (Bs.YZ);

B_Ha.cc = (Bs.ZZ) .* (Bs.ZZ);
B_Ha.cd = (Bs.ZZ) .* (Bs.XY);
B_Ha.ce = (Bs.ZZ) .* (Bs.XZ);
B_Ha.cf = (Bs.ZZ) .* (Bs.YZ);

B_Ha.dd = (Bs.XY) .* (Bs.XY);
B_Ha.de = (Bs.XY) .* (Bs.XZ);
B_Ha.df = (Bs.XY) .* (Bs.YZ);

B_Ha.ee = (Bs.XZ) .* (Bs.XZ);
B_Ha.ef = (Bs.XZ) .* (Bs.YZ);

B_Ha.ff = (Bs.YZ) .* (Bs.YZ);



end



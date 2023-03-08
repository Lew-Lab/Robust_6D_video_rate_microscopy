function locData = analyzeViaMVRestimator_v7(rawData,b,Bstruct,systemPar,RoSEreg,fName,varargin)
% b = background
% varargin - waitbar
s = opt2struct(varargin);

if isfield(s,'cost')
    if strcmp(s.cost,'Gauss')
        gaussCost = 1;
    else
        gaussCost = 0;
    end
else
    gaussCost = 0;
end

Blist = Bstruct.Blist;
BgradxList = Bstruct.BgradxList;
BgradyList = Bstruct.BgradyList;
BgradzList = Bstruct.BgradzList;
FPSF_r = Bstruct.FPSF_r;
FPSF_a = Bstruct.FPSF_a;
BsList = Bstruct.BsList;
B_HaList = Bstruct.B_HaList;
sumNormList = Bstruct.sumNormList;
zRoSE = systemPar.zRoSE;
zStep = systemPar.zStep;
zList = systemPar.zList;
PSFsz = systemPar.PSFsz;
pxSize = systemPar.pxSize;
NA = systemPar.NA;
offset = systemPar.offset;

locData = [];
obj.pixelUpsample = 1;
obj.pixelSize = pxSize;

if numel(b) == 8
    b = kron(b,ones(PSFsz));
end

for chNum = 1:8
    if size(b,1) == PSFsz % rotates background by k*90 degrees, where k is 5 - remainder of (channel number)//4
        b_ch{chNum} = rot90(b(:,PSFsz*(chNum-1)+(1:PSFsz)),5-rem(chNum,4));
    elseif numel(b) == 1
        b_ch{chNum} = b;
    end
end
totFrame = size(rawData,3);

D = parallel.pool.DataQueue;
if isfield(s, 'waitbar')
    wb = waitbar(0,s.waitbar,'estimate 0% complete');
    UpdateWaitbar(totFrame,wb);
    afterEach(D, @UpdateWaitbar);
end

for groupNum = 1:ceil(totFrame/500)
    frameList = ((groupNum-1)*500+1):min([groupNum*500,totFrame]);
    parfor frameNum = frameList
        img = double(rawData(:,:,frameNum))-offset;
        img(:,1:end/2) = img(:,1:end/2)*.49;
        img(:,end/2+1:end) = img(:,end/2+1:end)*.47;
        img(img<0) = 0;
        locDataFrame = [];
        for chNum = 1:8
            SMLM_img = rot90(img(:,PSFsz*(chNum-1)+(1:PSFsz)),5-rem(chNum,4));
            
            % Channels 1 - 4 handle radially (x) polarized light.
            % Channels 5 - 8 handle azimuthally (y) polarized light.
            if chNum <= 4
                FPSF = FPSF_r;
            else
                FPSF = FPSF_a;
            end
            
            % Temporary localization data obtained from 2D RoSE for each
            % channel number
            [~,~,locDataTmp1] = RoSE2D(obj, SMLM_img, b_ch{chNum}, FPSF, 'regVal', RoSEreg);
            if ~isempty(locDataTmp1)
                locDataTmp1(:,[2,3]) = locDataTmp1(:,[2,3])-pxSize;
                locDataTmp1(:,1) = frameNum;
                locDataTmp1(:,5) = chNum;
                tau = 280*mean(b_ch{chNum}(:))./locDataTmp1(:,4);
                locDataTmp1(:,6) = 1e-9*sqrt(400^2./locDataTmp1(:,4).*(1+4*tau+sqrt(2*tau./(1+4*tau))));
                locDataTmp2 = locDataTmp1;
                switch chNum
                    case 2
                        locDataTmp1(:,2) = locDataTmp2(:,3);
                        locDataTmp1(:,3) = -locDataTmp2(:,2);
                    case 3
                        locDataTmp1(:,2) = -locDataTmp2(:,2);
                        locDataTmp1(:,3) = -locDataTmp2(:,3);
                    case 4
                        locDataTmp1(:,2) = -locDataTmp2(:,3);
                        locDataTmp1(:,3) = locDataTmp2(:,2);
                    case 6
                        locDataTmp1(:,2) = locDataTmp2(:,3);
                        locDataTmp1(:,3) = -locDataTmp2(:,2);
                    case 7
                        locDataTmp1(:,2) = -locDataTmp2(:,2);
                        locDataTmp1(:,3) = -locDataTmp2(:,3);
                    case 8
                        locDataTmp1(:,2) = -locDataTmp2(:,3);
                        locDataTmp1(:,3) = locDataTmp2(:,2);
                end
            end
            locDataFrame = [locDataFrame;locDataTmp1];
        end
        if ~isempty(locDataFrame)
            locDataFrame(locDataFrame(:,4)<max([mean(b(:))*20,100]),:) = [];
        end
        if ~isempty(locDataFrame)
            groupCenterReg = 5;
            B_mle = 0;
            imgCrop = [];
            bCrop = [];
            for chNum = 1:8
                imgCrop = [imgCrop,img(6:end-5,(chNum-1)*PSFsz+(6:PSFsz-5))];
                if numel(b) ~= 1
                    bCrop = [bCrop,b(6:end-5,(chNum-1)*PSFsz+(6:PSFsz-5))];
                else
                    bCrop = b;
                end
            end
            img_mle = reshape(imgCrop,1,[]);
            b_mle = reshape(bCrop,1,[]);
            while min(abs(eig(B_mle*B_mle'))) < eps && groupCenterReg < 8
                x = []; y = []; z = [];
                mux = []; muy = []; muz = []; gamma = [];
                xRes = []; yRes = []; zRes = []; 
                xInd = []; yInd = []; zInd = [];
                cGroup = initXYZest(locDataFrame,frameNum,groupCenterReg,pxSize);
                B_mle = [];
                for locCount = 1:size(cGroup,1)
                    c = cGroup(locCount,:);
                    xInit = c(1);
                    yInit = c(2);
                    zInit = 1/(1.37 - NA/5)*c(3)+zRoSE(2); % approx equation
                    [B,xIndTmp,yIndTmp,zIndTmp] = shiftBasis(xInit,yInit,zInit,pxSize,Blist,BgradxList,BgradyList,BgradzList,zStep,zList);
                    xInd(locCount) = xIndTmp;
                    yInd(locCount) = yIndTmp;
                    zInd(locCount) = zIndTmp;
                    B_mle = [B_mle;B];
                end
                m = ((img_mle-b_mle)*B_mle')/(B_mle*B_mle');
                for iter = 1:10
                    B_mle = [];
                    for locCount = 1:size(cGroup,1)
                        mTmp = m((locCount-1)*15+(1:15));
                        [xTmp,yTmp,zTmp,xResTmp,yResTmp,zResTmp] = ...
                            refineLoc(mTmp,pxSize,zStep,zList,xInd(locCount),yInd(locCount),zInd(locCount));
                        x(locCount) = xTmp;
                        y(locCount) = yTmp;
                        z(locCount) = zTmp;
                        xRes(locCount) = xResTmp;
                        yRes(locCount) = yResTmp;
                        zRes(locCount) = zResTmp;
                        [B,xIndTmp,yIndTmp,zIndTmp] = shiftBasis(xTmp,yTmp,zTmp,...
                            pxSize,Blist,BgradxList,BgradyList,BgradzList,zStep,zList);
                        xInd(locCount) = xIndTmp;
                        yInd(locCount) = yIndTmp;
                        zInd(locCount) = zIndTmp;
                        B_mle = [B_mle;B];
                    end
                    if min(abs(eig(B_mle*B_mle'))) < eps
                        break;
                    end
                    if max([abs(xRes)/pxSize,abs(yRes)/pxSize,abs(zRes)/zStep]) < .5
                        break;
                    end
                    m = ((img_mle-b_mle)*B_mle')/(B_mle*B_mle');
                end
                groupCenterReg = groupCenterReg + 1;
            end
            if gaussCost
                if isfield(s,'gaussvar')
                    gaussVar = s.gaussvar;
                else
                    gaussVar = 25;
                end
                m = (mle_fista_Gauss(B_mle',img_mle',b_mle',m',gaussVar))';
            else
                m = (mle_fista(B_mle',img_mle',b_mle',m'))';
            end
            for locCount = 1:size(cGroup,1)
                mTmp = m((locCount-1)*15+(1:15));
                [xTmp,yTmp,zTmp] = refineLoc(mTmp,pxSize,zStep,zList,xInd(locCount),yInd(locCount),zInd(locCount));
                x(locCount) = xTmp;
                y(locCount) = yTmp;
                z(locCount) = zTmp;
                [muxTmp, muyTmp, muzTmp, gammaTmp] = secondM2SymmConeWeighted(BsList{zInd(locCount)},B_HaList{zInd(locCount)},...
                    sumNormList{zInd(locCount)},mTmp(1:6)/sum(mTmp(1:3)),sum(mTmp(1:3)),mean(b_mle(:)));
                mux(locCount) = muxTmp;
                muy(locCount) = muyTmp;
                muz(locCount) = muzTmp;
                gamma(locCount) = gammaTmp;
            end
            if ~isempty(x)
                locData = [locData,[frameNum*ones(1,length(x));x;y;z;mux;muy;muz;gamma;reshape(m,15,[])]];
            end
        end
        send(D,frameNum)
    end
    if ~isempty(fName)
        save([fName,'_locData.mat'],'locData');
    end
end
if ~isempty(locData)
    locData = locData(1:14,:);
end
locData(15,:) = sum(locData(9:11,:),1);
locData(:,locData(15,:)<max([mean(b(:))*50,200])) = []; % photon threshold

end

function [x,y,z,xRes,yRes,zRes] = refineLoc(m,pxSize,zStep,zList,xInd,yInd,zInd)
xRes = mean(m(7:9))./mean(m(1:3))*1e-7;
yRes = mean(m(10:12))./mean(m(1:3))*1e-7;
zRes = mean(m(13:15))./mean(m(1:3))*1e-7;

if abs(xRes)>pxSize*3
    xRes = sign(xRes)*pxSize*3;
end
if abs(yRes)>pxSize*3
    yRes = sign(yRes)*pxSize*3;
end
if abs(zRes)>zStep*3
    zRes = sign(zRes)*zStep*3;
end
x = xInd*pxSize + xRes;
y = yInd*pxSize + yRes;
z = (zInd-1)*zStep + zList(1) + zRes;
end

function UpdateWaitbar(totFrame,wb)
persistent TOTAL COUNT H
if nargin == 2
    TOTAL = totFrame;
    COUNT = 0;
    H = wb;
else
    COUNT = COUNT + 1;
    waitbar(COUNT/TOTAL, H, ['estimate ',num2str(round(COUNT/TOTAL*100)),'% complete']);
end
end
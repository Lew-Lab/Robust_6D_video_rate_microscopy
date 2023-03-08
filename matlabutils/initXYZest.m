function [cGroup,cGroupNoFilter,cFin] = initXYZest(loc,frameNum,tol,pxSize)

% cGroup - final group
% cGroupNoFilter - without sorting and deleting by score
% cFin - after matching two channels

for chNum = 1:8
    locTmp1 = loc(loc(:,5) == chNum,:);
    try
        locTmp1 = locTmp1(locTmp1(:,1)==frameNum,:);
        locTmp2{chNum} = locTmp1(:,[2,3,6]);
    catch
        continue;
    end
end
loc = RoSE2ChReg(locTmp2);

cFin = [];

for ch1 = 1:8
    for ch2 = (ch1+1):8
        c = findCenter2ch(loc{ch1},loc{ch2},ch1,ch2,tol,pxSize);
        cFin = [cFin;c];
    end
end

cGroup = [];
if ~isempty(cFin)
    [cGroup,cGroupNoFilter] = groupCenter(cFin,tol,pxSize);
end
for chNum = 1:8
    locTmp1 = loc{chNum};
    locTmp1 =[locTmp1,(1:size(locTmp1,1))'];
    if ~isempty(cFin)
        indTmp = cGroup(:,chNum+5);
        indTmp(isnan(indTmp)) = [];
        locTmp1(indTmp,:) = [];
    end
    spotInd = nan(size(locTmp1,1),8);
    spotInd(:,chNum) = locTmp1(:,4);
    cGroupTmp = [locTmp1(:,1),locTmp1(:,2),zeros(size(locTmp1,1),1),locTmp1(:,3),zeros(size(locTmp1,1),1),spotInd];
    cGroup = [cGroup;cGroupTmp];
end
cGroup = cGroup(:,1:3);

end

function [cGroup,cGroupNoFilter] = groupCenter(c,tol,pxSize)

cGroup = [];

c = [c,(1:size(c,1))'];

k = 1;
for i = 1:8
    for j = i+1:8
        indList(k,:) = [i,j];
        k = k + 1;
    end
end

for chPairNum1 = 1:28
    cTmp1 = c(c(:,5)==chPairNum1,:);
    for cTmp1Ind = 1:size(cTmp1,1)
        x = cTmp1(cTmp1Ind,1);
        y = cTmp1(cTmp1Ind,2);
        z = cTmp1(cTmp1Ind,3);
        prec = cTmp1(cTmp1Ind,4);
        sumx = x/prec;
        sumy = y/prec;
        sumz = z/prec;
        sumw = 1/prec;
        score = 1;
        for i = 6:13
            if ~isnan(cTmp1(cTmp1Ind,i))
                spotNum(i-5) = cTmp1(cTmp1Ind,i);
            else
                spotNum(i-5) = nan;
            end
        end
        c(c(:,14) == cTmp1(cTmp1Ind,14),:) = [];
        for chPairNum2 = chPairNum1+1:28
            cTmp2 = c(c(:,5)==chPairNum2,:);
            cTmp2Tmp = cTmp2;
            deleteCount = 0;
            for cTmp2Ind = 1:size(cTmp2Tmp,1)
                for i = 6:13
                    if ~isnan(cTmp2Tmp(cTmp2Ind,i)) && ~isnan(spotNum(i-5)) && cTmp2Tmp(cTmp2Ind,i)~=spotNum(i-5)
                        cTmp2(cTmp2Ind-deleteCount,:) = [];
                        deleteCount = deleteCount + 1;
                        break;
                    end
                end
            end
            [dist,cTmp2Ind] = min((cTmp2(:,1)-x).^2+(cTmp2(:,2)-y).^2+(cTmp2(:,3)-z).^2);
            dist = sqrt(dist);
            if dist < tol*(prec + cTmp2(cTmp2Ind,4)) + pxSize*2
                c(c(:,14) == cTmp2(cTmp2Ind,14),:) = [];
                sumx = sumx + cTmp2(cTmp2Ind,1)/cTmp2(cTmp2Ind,4);
                sumy = sumy + cTmp2(cTmp2Ind,2)/cTmp2(cTmp2Ind,4);
                sumz = sumz + cTmp2(cTmp2Ind,3)/cTmp2(cTmp2Ind,4);
                sumw = sumw + 1/cTmp2(cTmp2Ind,4);
                x = sumx/sumw;
                y = sumy/sumw;
                z = sumz/sumw;
                prec = prec*cTmp2(cTmp2Ind,4)/(prec+cTmp2(cTmp2Ind,4))*sqrt(2);
                score = score + 1;
                for i = 6:13
                    if ~isnan(cTmp2(cTmp2Ind,i))
                        spotNum(i-5) = cTmp2(cTmp2Ind,i);
                    end
                end
            end
        end
        cGroupTmp = [x,y,z,prec,score,spotNum];
        cGroup = [cGroup;cGroupTmp];
    end
end

cGroupNoFilter = cGroup;

% [~,ind] = sort(cGroup(:,5),'descend');
% cGroup = cGroup(ind,:);
cGroup = sortrows(cGroupNoFilter,[-5,4]);

for i = 1:size(cGroup,1)
    for spotInd = 6:13
        if ~isnan(cGroup(i,spotInd))
            cGroupTmp = cGroup(i+1:end,:);
            cGroupTmp(cGroupTmp(:,spotInd)==cGroup(i,spotInd),spotInd) = nan;
            cGroup(i+1:end,:) = cGroupTmp;
        end
    end
end

deleteCount = 0;
for i = 1:size(cGroup,1)
    if all(isnan(cGroup(i-deleteCount,6:13)))
        cGroup(i-deleteCount,:) = [];
        deleteCount = deleteCount+1;
    end
end

end



function c = findCenter2ch(posCh1,posCh2,ch1,ch2,tol,pxSize)

% posCh1: N1 x 3, each column represents x, y, loc prec
% posCh2: N2 x 3, each column represents x, y, loc prec

% c: N x 12 fitted center based on spot list from ch1 and ch2
% each column represents x, y, z, loc prec, spot index from list

x = [];
y = [];
z = [];
p = [];
chInd = [];
spotInd = [];

indAdd = 0;
for i = 1:7
    indAdd(i+1) = indAdd(i)+8-i;
end
chIndTmp = ch2-ch1+indAdd(ch1);

zmax = 2e-6;
normFac = [2*sind(22.5),sqrt(2),2*sind(67.5),2,2*sind(67.5),sqrt(2),2*sind(22.5)];

for indCh1 = 0:7
    x0(indCh1+1) = sind(45*indCh1);
    y0(indCh1+1) = cosd(45*indCh1);
end

v_para = [x0(ch2)-x0(ch1), y0(ch2)-y0(ch1)];
v_perp = [y0(ch1)-y0(ch2), x0(ch2)-x0(ch1)];
v_para = v_para/norm(v_para);
v_perp = v_perp/norm(v_perp);

for indCh1 = 1:size(posCh1,1)
    for indCh2 = 1:size(posCh2,1)
        v = posCh2(indCh2,1:2) - posCh1(indCh1,1:2);
        prec1 = posCh1(indCh1,3);
        prec2 = posCh2(indCh2,3);
        if abs(dot(v_perp,v)) < tol*(prec1 + prec2) + pxSize*2
            zTmp = dot(v_para,v) / normFac(rem(ch2-ch1+8,8));
            xTmp1 = posCh1(indCh1,1) - x0(ch1)*zTmp;
            yTmp1 = posCh1(indCh1,2) - y0(ch1)*zTmp;
            xTmp2 = posCh2(indCh2,1) - x0(ch2)*zTmp;
            yTmp2 = posCh2(indCh2,2) - y0(ch2)*zTmp;
            xTmp = (xTmp1/prec1 + xTmp2/prec2)/(1/prec1+1/prec2);
            yTmp = (yTmp1/prec1 + yTmp2/prec2)/(1/prec1+1/prec2);
            pTmp = prec1*prec2/(prec1+prec2)*sqrt(2);
            spotIndTmp = nan(1,8);
            spotIndTmp(ch1) = indCh1;
            spotIndTmp(ch2) = indCh2;
            if abs(zTmp) < zmax
                x = [x;xTmp];
                y = [y;yTmp];
                z = [z;zTmp];
                p = [p;pTmp];
                chInd = [chInd;chIndTmp];
                spotInd = [spotInd;spotIndTmp];
            end
        end
    end
end

c = [x,y,z,p,chInd,spotInd];

end
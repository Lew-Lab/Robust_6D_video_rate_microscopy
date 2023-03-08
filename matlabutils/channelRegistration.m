function locGroup = channelRegistration(loc,N,frameNum,thresholdL)

locGroup = {};
% N: max number of localizations in that frame

dist2 = nan(N,N,8,8);
dist1 = nan(N,N,8,8);

for chNum = 1:8
    locTmp = loc(loc(:,5) == chNum,:);
    try
        locTmp = locTmp(locTmp(:,1)==frameNum,:);
        xx{chNum} = locTmp(:,2);
        yy{chNum} = locTmp(:,3);
        II{chNum} = locTmp(:,4);
    catch
        continue;
    end
end

xx = RoSE2ChReg(xx);
yy = RoSE2ChReg(yy);
II = RoSE2ChReg(II);
    
for m = 1:7
    for n = (m+1):8
        for chNum = 1:numel(xx{m})
            for j = 1:numel(xx{n})
                diffVec = [xx{m}(chNum)-xx{n}(j),yy{m}(chNum)-yy{n}(j)];
                dist1(chNum,j,m,n) = norm(diffVec);
                dist2(chNum,j,m,n) = abs(dot(diffVec,[sind(22.5*(n-m)+45*(m-1)),cosd(22.5*(n-m)+45*(m-1))]));
            end
        end
    end
end

% thresholdL = 100e-9;
dist1(dist2>thresholdL) = nan;

normFac = [2*sind(22.5),sqrt(2),2*sind(67.5),2,2*sind(67.5),sqrt(2),2*sind(22.5)]/(2*sind(22.5));
locNum = 1;
  
for m = 1:7
   for n = m+1:8
        for chNum = 1:numel(xx{m})
            for j = 1:numel(xx{n})
                score = 0;
                idx = nan(1,8);
                if ~isnan(dist1(chNum,j,m,n))
                    for k = 1:8
                        if k == m || k == n
                            distDiff = nan;
                        else
                            distDiff = abs(dist1(j,:,n,k)/normFac(abs(k-n))-dist1(chNum,j,m,n)/normFac(n-m));
                        end
                        if nanmin(distDiff) < thresholdL
                            score = score + 1;
                            [~,idx(k)] = nanmin(distDiff);
                        end
                        idx(m) = chNum;
                        idx(n) = j;
                    end
                    if score >= 2
                        locGroup{locNum} = nan(3,8);
                        for ch = 1:8
                            if ~isnan(idx(ch))
                                locGroup{locNum}(1,ch) = xx{ch}(idx(ch));
                                locGroup{locNum}(2,ch) = yy{ch}(idx(ch));
                                locGroup{locNum}(3,ch) = II{ch}(idx(ch));
                                dist1(idx(ch),:,ch,:) = nan;
                                dist1(:,idx(ch),:,ch) = nan;
                            end
                        end
                        locNum = locNum + 1;
                    end
                end
            end
        end
    end
end

end

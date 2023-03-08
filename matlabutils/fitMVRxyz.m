function cPar = fitMVRxyz(locGroup)

sideLength = [2*sind(22.5),sqrt(2),2*sind(67.5),2,2*sind(67.5),sqrt(2),2*sind(22.5)];
sideLengthx = [0,-sqrt(.5),-1,-sqrt(.5),0,sqrt(.5),1,sqrt(.5)];
sideLengthy = [-1,-sqrt(.5),0,sqrt(.5),1,sqrt(.5),0,-sqrt(.5)];

x = locGroup(1,:);
y = locGroup(2,:);
w = sqrt(locGroup(3,:));
d = nan(7);
for i = 1:7
    for j = i+1:8
        d(i,j-1) = norm([x(i)-x(j),y(i)-y(j)]/sideLength(abs(i-j)));
    end
end
if nanmean(y([1:3]))>nanmean(y(5:7)) || nanmean(x(3:5))>nanmean(x([1,7,8]))
    rInit = nanmean(reshape(d,1,[]));
else
    rInit = -nanmean(reshape(d,1,[]));
end
xInit = nanmean(x+sideLengthx*rInit);
yInit = nanmean(y+sideLengthy*rInit);
cPar = fminsearch(@(c)circleCost(x,y,w,c),[rInit,xInit,yInit]);
% cPar = [rInit,xInit,yInit];
% figure(); hold();
% c = cPar;
% plot(x,y,'x');
% plot(c(2)+c(1)*sin(linspace(-pi,pi,100)),c(3)+c(1)*cos(linspace(-pi,pi,100)),'color',[0,0,1]);

end





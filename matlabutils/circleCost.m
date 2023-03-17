function cost = circleCost(x,y,w,c)

r = c(1);
cx = c(2);
cy = c(3);

t = [135,225,315,405,90,180,270,360];

cost = 0;
for i = 1:8
    if ~isempty(x(i)) && ~isnan(x(i))
        cost = cost + norm([x(i),y(i)] - [cx+r*cosd(t(i)),cy+r*sind(t(i))])*w(i);
    end
end

end

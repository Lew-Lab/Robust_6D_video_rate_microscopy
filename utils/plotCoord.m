function zz = plotCoord(data1)

minV = min(data1,[],"all");
maxV = max(data1,[],"all");
if minV==maxV
    minV = maxV-1;
end
caxis([minV,maxV]);

zz = max(data1,[],"all")*40;
if zz<=1
    zz=1;
end
zlim([0,zz]);

plot3(sin(linspace(-pi,pi,100)),cos(linspace(-pi,pi,100)),zz*ones(1,100),'k','linewidth',.5);
plot3(sind(30)*sin(linspace(-pi,pi,100)),sind(30)*cos(linspace(-pi,pi,100)),zz*ones(1,100),'k','linewidth',.5);
plot3(sind(60)*sin(linspace(-pi,pi,100)),sind(60)*cos(linspace(-pi,pi,100)),zz*ones(1,100),'k','linewidth',.5);
for p = 0:30:150
    plot3(sind(p)*[-1,1],cosd(p)*[-1,1],zz*[1,1],'k','linewidth',.2);
end
set(gca,'visible','off');

% text(1.1*sind(0),1.1*cosd(0),zz,'$90^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle','rotation',0)
text(1.05*sind(60),1.05*cosd(60),zz,'$30^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle','rotation',-60)
text(1.05*sind(90),1.05*cosd(90),zz,'$\phi=0^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle','rotation',-90)
text(1.05*sind(120),1.05*cosd(120),zz,'$330^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle','rotation',-120)
text(1.05*sind(240),1.05*cosd(240),zz,'$210^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle','rotation',-240)
text(1.05*sind(300),1.05*cosd(300),zz,'$150^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle','rotation',-300)
text(sind(175)*sind(5),cosd(175)*sind(5),zz,'$0^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','left','VerticalAlignment','middle')
text(sind(175)*sind(30),cosd(175)*sind(30),zz,'$30^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','left','VerticalAlignment','middle')
text(sind(175)*sind(60),cosd(175)*sind(60),zz,'$60^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','left','VerticalAlignment','middle')
text(sind(175),cosd(175)*1.1,zz,'$\theta=90^\circ$','interpreter','latex','fontsize',10,'HorizontalAlignment','left','VerticalAlignment','middle')


end
clear; close all;

% load cMapPhase.mat;
% load 2000_5_z2_locData.mat; % best
% c = [210,-40,0];
% r = 350;
% Nz = 500;
% Nxy = 2500;
% locData = driftCorrection(locData,c,r,Nz,Nxy,200,200,3000,'off');
% locData = locData(:,abs(locData(2,:))<2000);
% locData = locData(:,abs(locData(3,:))<2000);
% locData_z2 = locData(:,locData(4,:)>=1000);

% load 2000_5_z3_locData.mat; % best
% locData(15,:) = sum(locData(9:11,:));
% locData(16:end,:) = [];
% c = [300,-50,0];
% r = 350;
% Nz = 500;
% Nxy = 2500;
% locData = driftCorrection(locData,c,r,Nz,Nxy,200,200,3000,'off');
% locData = locData(:,abs(locData(2,:))<2000);
% locData = locData(:,abs(locData(3,:))<2000);
% locData_z3 = locData(:,locData(4,:)>=1000);

load 2000_5_z1_locData.mat; % best
c = [420,200,0];
r = 350;
Nz = 500;
Nxy = 2500;

% locData = driftCorrection(locData,c,r,Nz,Nxy,zThreshold,xyThreshold,ROIsz,vis)
locData = driftCorrection(locData,c,r,Nz,Nxy,200,200,3000,'off');
locData = locData(:,abs(locData(2,:))<2000);
locData = locData(:,abs(locData(3,:))<2000);
locData_z1 = locData(:,locData(4,:)<=1000);

locData = [locData_z1,locData_z3];

locData = locData(:,locData(15,:)>1500);

% figure();
% hold();
% fill(2000*[-1,-1,1,1,-1],2000*[-1,1,1,-1,-1]+400,'k');
% plotTheta(locData(2,:),locData(4,:),locData(5,:),locData(7,:),20,parula(256))
% axis image;
% ylim([-100,2400]);
% 
% figure();
% for i = 1:9
%     subplot(3,3,i); hold();
%     fill(2000*[-1,-1,1,1,-1],2000*[-1,1,1,-1,-1],'k');
%     idx = logical((locData(4,:)<i*100).*(locData(4,:)>=(i-1)*100));
%     plotPhi(locData(2,idx),locData(3,idx),locData(5,idx),locData(6,idx),20,cMapPhase);
%     colormap(cMapPhase);
%     axis image;
% end
% % 
% % figure();
% % hold();
% % fill(500*[-1,-1,1,1,-1],500*[-1,1,1,-1,-1]+400,'k');
% % scatter(locData(2,:),locData(4,:),[],acosd(locData(7,:)),'filled','SizeData',5,'MarkerFaceAlpha',.5);
% % axis image;
% % ylim([-100,900]);
% 
% omega = (3-sqrt(1+8*locData(8,:)))*pi;
% % figure();
% % histogram(omega)
% 
% % figure();
% 
% 
% 
% function plotPhi(x,y,mux,muy,r,cMap)
% p = plot([x+r*mux;x-r*mux],[y+r*muy;y-r*muy],'linewidth',1);
% phi = atand(muy./mux);
% cList = cell(length(x),1);
% for i = 1:length(phi)
%     cList{i} = cMap(1+round((length(cMap)-1)*(phi(i)+90)/180),:);
% end
% set(p,{'color'},cList);
% end
% 
% function plotTheta(x,z,mux,muz,r,cMap)
% p = plot([x+r*mux;x-r*mux],[z+r*muz;z-r*muz],'linewidth',1);
% theta = acosd(abs(muz));
% cList = cell(length(x),1);
% for i = 1:length(theta)
%     cList{i} = cMap(1+round((length(cMap)-1)*(theta(i))/90),:);
% end
% set(p,{'color'},cList);
% end
% 

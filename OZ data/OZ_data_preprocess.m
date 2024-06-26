clear; clc;

fname =  '2000_5_z1';
offset = load([fname '_offset.mat']).offset;
rawData = load([fname '_rawData.mat'],'-mat').rawData;
bg = load([fname '_bg.mat'],'-mat').b;
numFrame = size(rawData,3);
data = zeros(numFrame,8,size(rawData,1),size(rawData,1));

for i = 1:numFrame
    dataTmp = double(rawData(:,:,i)) - offset;
    % camera brightness correction
    dataTmp(:,1:end/2) = dataTmp(:,1:end/2)*.49;
    dataTmp(:,end/2+1:end) = dataTmp(:,end/2+1:end)*.47;

    temp = zeros(8,size(dataTmp,1),size(dataTmp,1));
    for j = 1:8 
        temp(j,:,:) = dataTmp(:,size(dataTmp,1)*(j-1)+1:size(dataTmp,1)*j)';
        data(i,j,:,:) = temp(j,:,:);
    end
    data_cat(:,:,i) = cat(2,squeeze(temp(1,:,:)),squeeze(temp(2,:,:)),...
        squeeze(temp(3,:,:)),squeeze(temp(4,:,:)),...
        squeeze(temp(5,:,:)),squeeze(temp(6,:,:)),...
        squeeze(temp(7,:,:)),squeeze(temp(8,:,:)));
end
%% 

dataSum = sum(data_cat,3);

image = zeros(8,size(dataSum,1),size(dataSum,1));
for i = 1:8 
    image(i,:,:) = dataSum(:,size(dataSum,1)*(i-1)+1:size(dataSum,1)*i)-bg(i);
end

%%

for i = 1:8
    image(i,:,:) = image(i,:,:)-bg(i);
end
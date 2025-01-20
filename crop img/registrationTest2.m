clear

fileFolder = 'C:\Users\cheng\OneDrive - Washington University in St. Louis\Desktop\Lew Lab\Robust-6D-video-rate-microscopy\experiment data\';
fNameRead = 'fluorescent beads_crop';
fNameFull = [fileFolder,fNameRead,'.tif'];
fNameSave = 'registration_check';

% load the csv file and calculate the relative position
locTmp = csvread(['registration check','.csv'],1,0);
locTmp = locTmp(:,2:3)/66.857;

color = ['r' 'g' 'b' 'g' 'y' 'k' 'c' 'm' 'r'];
figure
for i = 1:8
    con = locTmp(:,1)<141*i & locTmp(:,1)>141*(i-1);
    loc = locTmp(con,:);
    scatter(loc(:,1)-141*(i-1),loc(:,2),color(i))
    hold on
end

title('localization of fluorescent beads')
xlabel('x');
ylabel('y');
xlim([1 141])
ylim([1 141])

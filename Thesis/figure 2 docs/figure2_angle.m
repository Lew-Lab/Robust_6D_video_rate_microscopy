theta = [5:10:85];
theta_error_mean = zeros(9,1);
theta_error_std = zeros(9,1);
theta_error_25per = zeros(9,1);
theta_error_50per = zeros(9,1);
theta_error_75per = zeros(9,1);

for i = 1:9
    one = readmatrix(['theta ' num2str((i-1)*2) '.csv']);
    two = readmatrix(['theta ' num2str((i-1)*2+1) '.csv']);
    temp = [one(2:end,2);two(2:end,2)];
    theta_error_mean(i) = mean(temp);
    theta_error_std(i) = std(temp);
    theta_error_25per(i) = prctile(temp,25);
    theta_error_50per(i) = prctile(temp,50);
    theta_error_75per(i) = prctile(temp,75);
end

figure;
plot(theta,theta_error_50per);
hold on
plot(theta,theta_error_25per);
hold on
plot(theta,theta_error_75per);
hold on

%%
phi = [-170:20:170];
phi_error_mean = zeros(18,1);
phi_error_std = zeros(18,1);
phi_error_25per = zeros(18,1);
phi_error_50per = zeros(18,1);
phi_error_75per = zeros(18,1);

for i = 1:18
    one = readmatrix(['phi ' num2str((i-1)*2) '.csv']);
    two = readmatrix(['phi ' num2str((i-1)*2+1) '.csv']);
    temp = [one(2:end,2);two(2:end,2)];
    phi_error_mean(i) = mean(temp);
    phi_error_std(i) = std(temp);
    phi_error_25per(i) = prctile(temp,25);
    phi_error_50per(i) = prctile(temp,50);
    phi_error_75per(i) = prctile(temp,75);
end

figure;
plot(phi,phi_error_50per);
hold on
plot(phi,phi_error_25per);
hold on
plot(phi,phi_error_75per);
hold on

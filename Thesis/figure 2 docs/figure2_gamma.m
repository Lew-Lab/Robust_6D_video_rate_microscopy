gamma_bin = linspace(0,2*pi,31);
gamma_0 = readmatrix('gamma 0 est.csv');
gamma_0 = gamma_0(2:end,2);
gamma_0 = 4*pi^2*(3/(4*pi)-sqrt(9/(16*pi^2)-(1-gamma_0)/(2*pi^2)));
gamma_1 = readmatrix('gamma 1 est.csv');
gamma_1= gamma_1(2:end,2);
gamma_1 = 4*pi^2*(3/(4*pi)-sqrt(9/(16*pi^2)-(1-gamma_1)/(2*pi^2)));
gamma_half = readmatrix('gamma half est.csv');
gamma_half= gamma_half(2:end,2);
gamma_half = 4*pi^2*(3/(4*pi)-sqrt(9/(16*pi^2)-(1-gamma_half)/(2*pi^2)));

histogram(gamma_0,gamma_bin,Normalization="probability")
hold on
histogram(gamma_1,gamma_bin,Normalization="probability")
hold on
histogram(gamma_half,gamma_bin,Normalization="probability")
hold on
xline(mean(gamma_1))
hold on 
xline(mean(gamma_0))
hold on
xline(mean(gamma_half))

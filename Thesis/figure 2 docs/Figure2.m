zlayer = [1,2,3,4,5,6,7,8,9,10];

fp = [0.000656932
0.009898773
0.013467109
0.014557019
0.005098689
0.00296366
0.004912061
0.01843889
0.017558004
0.015512556
];

fn = [0
0.000238884
0.000156768
0.000515095
0.000544955
0.003023381
0.001463167
0.000373257
0
0
];

w1 = .75; 
barh(zlayer,fp,w1)
w2 = .5;
hold on
barh(zlayer,fn,w2)
hold off

legend('fp rate', 'fn rate')

%% gamma vs. gamma error

gamma = [0, .5, 1];

gamma_error_50per = [0.079983715	0.503098314	0.997203187];
gamma_error_25per = [0.054126489	0.437949321	0.921209225];
gamma_error_75per = [0.155152048	0.687969902	1];

plot(gamma,gamma_error_75per,'b',LineWidth=1);
hold on
plot(gamma,gamma_error_25per,'b',LineWidth=1);
hold on
plot(gamma,gamma_error_50per,'b',LineWidth=1);

%% gamma vs. brt error










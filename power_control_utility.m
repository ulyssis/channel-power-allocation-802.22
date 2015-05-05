%%%% the utility for power allocation
%%%% 
%%%% E: unit is Joule
%%%% R: transmission rate: unit is bits/sec
%%%% L: length of packet

E= 1;
R=1;
L=3;
W=6;
interference = 1;
h=1;
f = 0.5;

p= 0.01:0.01:10;
gamma = W/R * h*p/f; % sir
figure;
u = E*R./p.*(1-exp(-0.5*gamma)).^L;


plot(p, u);





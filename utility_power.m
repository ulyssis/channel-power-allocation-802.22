% W: spectrum bandwidth
% R: rate, bits/sec
% f: interferece 
% h: path gain
% E: engerny
% L: length of packet
%


W = 6e+6;
R = 1e+6;
% f = 1e-6;
% h = 1e-6;
% f = 3.0412e-09;
% h= 2.0408e-08;

E=1;
L=80;
t = 1;

p=0.1:0.1:10;
gamma = W*h*p/(R*f);
% u = E*R./p .* (1-exp(-0.5.*gamma)).^L-t*p;
u = E*R./p .* (1-exp(-0.5.*gamma)).^L;
plot(p, u);


[x,y] = meshgrid(0:1:100);
z = (x.^2+y.^2)./(x.*y);

mesh(x,y, z);
grid on;
axis on;
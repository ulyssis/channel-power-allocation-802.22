%to test the function
clf
axis equal
hold on
x1=1;
y1=1;
rc=1;
[x,y,z] = cylinder(rc,200);
plot(x(1,:) + x1, y(1,:) + y1, 'r');
for t=1:20 %loop until doing 1000 points inside the circle
[x y]=cirrdnPJ(x1,y1,rc);
plot(x,y,'x');
%pause(0.01) %if you want to see the point being drawn
end
% % round distribution
% function [x y] = cirrdnPJ(x1,y1,rc)
% %the function, must be on a folder in matlab path
% % randNum = 0.1+0.7*rand;
% % randNum = rand;
% % a = 2*pi*randNum;
% a = 2*pi*rand;
% % r = sqrt(rand);
% x = (rc*(0.1+rand*0.9))*cos(a)+x1;
% y = (rc*(0.1+rand*0.9))*sin(a)+y1;
% 
% % x = rc*cos(a)+x1;
% % y = rc*sin(a)+y1;
% end


% square distribution
function [x, y] = cirrdnPJ(x1,y1,rc)

a = rc*(2*rand-1);
x = x1 + a;
b = rc*(2*rand-1);
y = y1 + b;

% x = rc*cos(a)+x1;
% y = rc*sin(a)+y1;
end

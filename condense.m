function [condenseB] = condense(B, w)

n = size(B, 1)/w;
c = size(B, 2);
condenseB = zeros(n, c);

for i = 1 : w
    condenseB = condenseB + B((i-1)*n+1: i*n, 1:c);
end

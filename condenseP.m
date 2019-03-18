function [condenseP_CVX] = condenseP(P_CVX, n, c)

condenseP_CVX = zeros(n, c);

for i = 1: n

       condenseP_CVX(i, :) =  sum(P_CVX((i-1)*c+1: i*c, :), 1);
end

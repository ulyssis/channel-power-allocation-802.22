% Randomly assign power to SUs

function [P] = maximalPowerGive(n, c)
    
P = zeros(n*c ,c);
    basePower = 4;
    for su = 1:n
        for i = 1:c 
            P((su-1)*c + i, i) = basePower + floor(rand*6);
        end
    end
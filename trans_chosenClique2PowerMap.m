function [powerMap] =  trans_chosenClique2PowerMap(chosenChannel, chosenClique, n, c, P_CVX)

powerMap = zeros(n, c);
for i = 1:n
    if(chosenClique(i))
       powerMap(i, :) = P_CVX((i-1)*c + chosenChannel, :);
    end
end

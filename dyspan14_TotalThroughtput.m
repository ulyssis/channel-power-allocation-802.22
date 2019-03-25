%- just a template
% [B, updateFlag] = update(seq(i), w, B, P_CVX, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta);

% chosenClique is a  1Xn vector, elements are 0/1, 1 means the correspoding WBS is in the clique. 
% n: number of physical WBSs
% P_CVX is n X c
% Gtilde is n X c
% chosenClique is 1 X n
function [totalThroughput] = dyspan14_TotalThroughtput(chosenChannel, chosenClique, n, c, P_CVX, Gtilde, TVpower, delta, eta)


        %B(su, :) = P_CVX((realIndex-1)*c + chosenChannel, :);
        F = (chosenClique' * chosenClique ~= 0);
        F = F - eye(n);
        powerMap = trans_chosenClique2PowerMap(chosenChannel, chosenClique, n, c, P_CVX);
        InterferenceonAll = Gtilde.* F * sum(powerMap, 2)*eta + delta;
        part1 = log2(1+ sum(powerMap, 2)./InterferenceonAll);

        totalThroughput = sum(part1);

        

%- just a template
% [B, updateFlag] = update(seq(i), w, B, P_CVX, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta);

% chosenClique is a  1Xn vector, elements are 0/1, 1 means the correspoding WBS is in the clique. 
function [totalThroughput] = dyspan14_TotalThroughtput(chosenChannel, chosenClique, P_CVX, Gtilde, TVpower, delta, SUcellRadius, pathlossfactor, eta)

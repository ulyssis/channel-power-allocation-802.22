function [utility] = U(i, channel, B, P, Gtilde, m, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta)

n = size(B, 1);
c = size(B, 2);
B(i, :) = P((i-1)*c + channel, :);

F = (B* B' ~= 0);
F = F - eye(size(B, 1));
    InterferenceonAll = Gtilde.* F * sum(B, 2)*eta + sum(GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) * TVpower, 2) + delta;
    utility = InterferenceonAll(i)/(sum(B(i,:))/SUcellRadius^(-pathlossfactor));
    

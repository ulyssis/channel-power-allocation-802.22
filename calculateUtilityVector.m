% Obtain a n x 1 utility function 
function [tem] = calculateUtilityVector(B, n, Gtilde, delta)
F = (B* B' ~= 0);   % F illustrates the interferce relations
F = F - eye(n);
    % The initial utility
    InterferenceonAll = Gtilde.* F * sum(B, 2) + delta;
    part1 = InterferenceonAll./ sum(B, 2);   % n X 1
    
    dummysumB = ones(n,1)*sum(B, 2)';    % duplicate B into a nxn matrix
    part2 = (Gtilde.*F.*dummysumB)./(sum(B, 2)*ones(1,n));  % n X n
    tem = part1 + sum(part2)';

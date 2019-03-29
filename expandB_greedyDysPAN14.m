function [expandedB] = expandB_greedyDysPAN14(B)

    n = size(B, 1);
    c = size(B, 2);

    tol = 1.e-6;
    B(B<tol) = 0;

    expandedB = zeros(n*c, c);
    
    for i = 1:n
        oneRow = B(i, :);
        oneBlock = diag(oneRow);
        expandedB ((i-1)*c+1:i*c, :) = oneBlock;
    end



% if(size(expandB, 1) ~= n*w)
%    stop =1; 
% end
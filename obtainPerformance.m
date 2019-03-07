% sum up the utilites over all nodes
% Obtain a n x 1 utility function 

% Obtain the following items:
% sum of utilities over all n users.
% Average interference
% Average transmssion power
% average SINR

function [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(w, B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu)
    F = (B* B' ~= 0);   % F illustrates the interferce relations
    F = F - eye(n);

        interfTV = sum(GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) *TVpower, 2);  % interference on every users

        InterferenceonAll = Gtilde.* F * repmat(Gtilde.* F, w, w) * sum(B, 2) + interfTV + delta;
    
    sumUtility = sum(InterferenceonAll./ ((sum(B, 2))*SUcellRadius^(-pathlossfactor)));    % f/p
    
    averageI = mean(InterferenceonAll); % average interference
    averageP = mean(sum(B, 2));
    averageSINR = mean(10*log10((sum(B, 2))*SUcellRadius^(-pathlossfactor)./InterferenceonAll));
%     averageSINR = mean(sum(B, 2)*SUcellRadius^(-pathlossfactor)./InterferenceonAll); % numerical, without log10
    stdSINR = std(10*log10((sum(B, 2))*SUcellRadius^(-pathlossfactor)./InterferenceonAll));
%     stdSINR = std(10*log10(sum(B, 2)*SUcellRadius^(-pathlossfactor)./InterferenceonAll));
% disp(sum(B, 2)./(InterferenceonAll+ interfTV + delta));


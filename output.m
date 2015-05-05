% print out the snrRatio vector which is n x 1

function     snrRatio = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor, eta)
F = (B* B' ~= 0);   % F illustrates the interferce relations
    F = F - eye(n);
    % The initial utility
    InterferenceonAll = Gtilde.* F * sum(B, 2) + sum(GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) *TVpower, 2) + delta;
%     snrRatio = 10*log10(sum(B, 2)*SUcellRadius^(-pathlossfactor)./InterferenceonAll);
    snrRatio = 10*log10(sum(B, 2)*SUcellRadius^(-pathlossfactor)./InterferenceonAll);
     disp(snrRatio);

    
    
    
    
    %     disp('Capacity:');
    %     disp(log10(sum(B, 2)./InterferenceonAll));
function [resultedInterference, exceedInterferenceBound] = checkResultedInference(B_scheme2Centralized, n, m, GtildeAll, infBound)



        
        % check the sum of interference, which include the  to the TV contour which works on channel
        % i.
        InfToTV = GtildeAll(n+1:n+m, 1:n)'.* B_scheme2Centralized;
        resultedInterference = sum(InfToTV);
        exceedInterferenceBound = sum(InfToTV) > infBound;


            
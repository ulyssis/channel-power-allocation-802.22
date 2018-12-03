   %% random channel allocation
   function [initialB] = generateInitialChannelAllocationB_for_schemeI(n,c, P)

        % Initialize channels asignment randomly
        B = zeros(n, c);
        doagain=1;
            while (doagain)
                for i = 1 : n
                   B(i, :) = P((i-1)*c + floor((1+c * rand)), :);        
                end
                doagain=0;
                for i=1:c
                   if (nnz(B(:,i))==1)
                       doagain=1;
                       break;
                   end
                end
            end

       initialB = B;       % record the initial B for selfishUpdate
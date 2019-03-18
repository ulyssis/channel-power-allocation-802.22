function [expandB] = expandB(B, w)

    n = size(B, 1);

    c = size(B, 2);
    
if(w ~=1)

    expandB = zeros(n*w, c);
    tol = 1.e-6;
    
    for i = 1:n

        Vindex = find(B(i, :) > tol);

        for j = 1: size(Vindex, 2)
            row = zeros(1, c);
            row(Vindex(j)) = 1;
            expandB((j-1)*n + i, :) = B(i, :).* row;	
        end
    end
else
    expandB = B;
end

if(size(expandB, 1) ~= n*w)
   stop =1; 
end
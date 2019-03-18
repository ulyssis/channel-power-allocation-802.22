function [realIndex] = getRealIndexFromTheExpanded(su, n, w)

for i = 1: w
   if(su > (i-1)*n/w && su <= i*n/w)
        realIndex = su - (i-1)*n/w;
   end
end


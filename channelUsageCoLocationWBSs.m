
% when a WBS looks for a channel, it should avoid using the channels being
% used the co-location WBSs.

function [thisChannelIsNotUsedByCOLocationWBSs] = channelUsageCoLocationWBSs(k, n, w, su, B)

c = size(B, 2);

for i = 1: w
   if(su > (i-1)*n && su <= i*n)
        realIndex = su - (i-1)*n;
   end
end


channelUsageCoLocationWBSs = zeros(w-1, c);
temp = 1;
for i = realIndex : n/w : n
    if(i~=su)
        channelUsageCoLocationWBSs(temp, :) = B(i, :);
        temp = temp + 1;
    end
end


if(nnz(channelUsageCoLocationWBSs(:, k)) == 0)
    thisChannelIsNotUsedByCOLocationWBSs = 1;
else
    thisChannelIsNotUsedByCOLocationWBSs = 0;
end
%output: a 1 X n vector, 0/1

function [cliqueUsingTheChannel, cliqueHavingTheChannel] = dyspan14_findCliqueUsingTheChannel(discussedChannel, n, channelAllocation, availableChannelsAllWBSs)

cliqueUsingTheChannel = zeros(1, n);
cliqueHavingTheChannel = zeros(1, n);

for i = 1: n
    
    thisWbsisUsing = ismember(discussedChannel, find(channelAllocation(i, :)));
    if(thisWbsisUsing)
        cliqueUsingTheChannel(i) = 1;
    end
    
    
    thisWbsHaving = ismember(discussedChannel, find(cell2mat(availableChannelsAllWBSs(i))));
    if(thisWbsHaving)
        cliqueHavingTheChannel(i) = 1;
    end
    
end
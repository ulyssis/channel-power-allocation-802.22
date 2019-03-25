%- just a template
function [channelAllocation] = dyspan14_GreedyAssign(n, c, P_CVX, Gtilde, channelAllocation, TVpower, delta, eta)


% 1: n cell, each element is a vector 1: c
availableChannelsAllWBSs = dyspan14_createReservedChannelsAllWBSs(n, c);

unionAvailableChannels = cell2mat(availableChannelsAllWBSs);

while(nnz(unionAvailableChannels))
    
    % Greedy algorithm
    
    % iterate all WBSs
    for i = 1:n
        
        % get the best channel: Cmax
        Cmax = dyspan14_BestChannel(i, n, c, P_CVX, Gtilde, channelAllocation, availableChannelsAllWBSs, TVpower, delta, eta);
        
        % assign the best channel to the current WBS
        channelAllocation(i, Cmax) = 1;
        
        % remove the best channel from the available channels.
        availableChannelsOnWBSi = cell2mat(availableChannelsAllWBSs(i));
        availableChannelsOnWBSi(Cmax) = 0;
        availableChannelsAllWBSs(i) = {availableChannelsOnWBSi};

        % remove the best channel from neighbors' available channels
        availableChannelsAllWBSs = dyspan14_RemoveNeighbors(n, c, Cmax, channelAllocation, availableChannelsAllWBSs, P_CVX, Gtilde, TVpower, delta, eta);
    end

    % update A, to check whether it is empty
    unionAvailableChannels = cell2mat(availableChannelsAllWBSs);
    
   
end


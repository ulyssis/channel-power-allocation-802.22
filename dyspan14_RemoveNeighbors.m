function  [availableChannelsAllWBSs] = dyspan14_RemoveNeighbors(n, c, Cmax, channelAllocation, availableChannelsAllWBSs, P_CVX, Gtilde, TVpower, delta, eta, SUcellRadius, pathlossfactor)

[cliqueUsingTheChannel, cliqueHavingTheChannel] = dyspan14_findCliqueUsingTheChannel(Cmax, n, channelAllocation, availableChannelsAllWBSs);


D1 = cliqueUsingTheChannel;
D2 = cliqueHavingTheChannel;

WbsIndiceHavingTheChannel = find(D2);

for elm = WbsIndiceHavingTheChannel
    
    tem = D1;
    tem(elm) =1;
   TT1 = dyspan14_TotalThroughtput(Cmax, tem, n, c, P_CVX, Gtilde, TVpower, delta, eta, SUcellRadius, pathlossfactor);
   TT2 = dyspan14_TotalThroughtput(Cmax, D1, n, c, P_CVX, Gtilde, TVpower, delta, eta, SUcellRadius, pathlossfactor);
   
   if(TT1 < TT2)     
       availableChannelOnElm = availableChannelsAllWBSs{elm};
       availableChannelOnElm(Cmax) = 0;
       availableChannelsAllWBSs(elm) = {availableChannelOnElm};
   end
end
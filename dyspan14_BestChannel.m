%- just a template
% Obtain the best channel.
% according to the pseudo code in ALg2 in the paper.
function [Cmax] = dyspan14_BestChannel(i, n, c, P_CVX, Gtilde, channelAllocation, availableChannelsAllWBSs, TVpower, delta, eta, SUcellRadius, pathlossfactor)

f = - 1000;

availableChannelsOnWBSi = availableChannelsAllWBSs{i};

ChannelIndice = find(availableChannelsOnWBSi);

if(isempty(ChannelIndice))
 
   return;
else
    for channelInTest = ChannelIndice

        [cliqueUsingTheChannel, cliqueHavingTheChannel] = dyspan14_findCliqueUsingTheChannel(channelInTest, n, channelAllocation, availableChannelsAllWBSs);

        D1 = cliqueUsingTheChannel;
        D1(i) =1;

        D2 = cliqueHavingTheChannel;    
        WbsIndiceHavingTheChannel = find(D2);

        for elm = WbsIndiceHavingTheChannel
            tem = D1;
            tem(elm) =1;
            TT1 = dyspan14_TotalThroughtput(channelInTest, tem, n, c, P_CVX, Gtilde, TVpower, delta, eta, SUcellRadius, pathlossfactor);
            TT2 = dyspan14_TotalThroughtput(channelInTest, D1, n, c, P_CVX, Gtilde, TVpower, delta, eta, SUcellRadius, pathlossfactor);
            if(TT1 > TT2)
                D1(elm) = 1;
            end
        end

        TT  = dyspan14_TotalThroughtput(channelInTest, D1, n, c, P_CVX, Gtilde, TVpower, delta, eta, SUcellRadius, pathlossfactor);
        if(TT > f)
            f = TT;
            Cmax = channelInTest;
        end

    end

end
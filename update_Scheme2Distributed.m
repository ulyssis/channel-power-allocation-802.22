% update_Scheme2Distributed, update channel

% su: index of current secondary user; 

function [B, updateFlag] = update_Scheme2Distributed(su, B, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta, infBound, POperation)
workingChannelProfile = B(su,:); % current channel usage

n = size(B, 1);
c = size(B, 2);   %   numeracally equal to c, 

UtilitySU = zeros(1, c);

atLeastAChannelIsAllowedToUse = 0;
for i = 1: c
        channelUsage = zeros(1, c);
        channelUsage(i) = POperation;
        B(su, :) = channelUsage;
        
        % check the sum of interference, which include the  to the TV contour which works on channel
        % i.
        InfToTV = GtildeAll(n+1:n+m, 1:n)'.* (B~=0) * POperation;
        if(sum(InfToTV(:,i)) > infBound)
            UtilitySU(i) = -1;
            B(su, i) = 0;
        else
            atLeastAChannelIsAllowedToUse = 1;

            % calculate the utility w.r.t. channel
            F = (B* B' ~= 0);
            F = F - eye(size(B, 1));
            InterferenceonAll = Gtilde.* F * sum(B, 2)*POperation*eta + delta/(n/c);
            part1 = InterferenceonAll(su)/(sum(B(su,:))*SUcellRadius^(-pathlossfactor));
            UtilitySU(i) = part1;
        end

end


if( ~any(workingChannelProfile))    % if true: WBS was idle before this update
    
    if(~atLeastAChannelIsAllowedToUse)
        updateFlag = 0;
        B(su, :) = zeros(1, c);
        return
    else
        updateFlag = 1;
        minUtility = min(UtilitySU(UtilitySU>0)); % minimum positive utility
        channel = find(UtilitySU == minUtility, 1);
        B(su, :) = zeros(1, c);
        B(su, channel) = POperation;
    end
else
    minUtility = min(UtilitySU(UtilitySU>0));   % minimum positive utility
    tem2 = find(UtilitySU == minUtility, 1);                                % obtain the channel which achieves the smallest utility
    workingChannel = find(workingChannelProfile, 1);
    if (tem2 ~= workingChannel && UtilitySU(workingChannel) > minUtility)   % There is another channel leading to smaller utility
        B(su, :) = zeros(1, c);
        B(su, tem2) = POperation;
        updateFlag = 1;
    else
        B(su, :) = zeros(1, c);
        B(su, workingChannel) = POperation;
        updateFlag = 0;
    end    
end








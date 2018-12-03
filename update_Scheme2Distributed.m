% update_Scheme2Distributed, update channel

% su: index of current secondary user; 

function [B, updateFlag] = update_Scheme2Distributed(su, B, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta, infBound, POperation)

workingChannel = find(B(su,:)); % which channel is currently used? or not at all?

n = size(B, 1);
c = size(B, 2);   %   numeracally equal to c, 

UtilitySU = zeros(1, c);

for i = 1: c
        channelUsage = zeros(1, c);
        channelUsage(i) = POperation;
        B(su, :) = channelUsage;
        
        % check the interference to the TV
        InfToTV = (GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) * TVpower)./((B+(B==0))*SUcellRadius^(-pathlossfactor));
        if(sum(InfToTV(:,i)) > Inf)
            UtilitySU(i) = 0;
           break; 
        end
        
        
        F = (B* B' ~= 0);
        F = F - eye(size(B, 1));
        InterferenceonAll = Gtilde.* F * sum(B, 2)*POperation*eta + delta/(n/c);
        part1 = InterferenceonAll(su)/(sum(B(su,:))*SUcellRadius^(-pathlossfactor));

        
        UtilitySU(i) = part1;
end


if( isempty(workingChannel))    % WBS switchs from idle to xxxx state
    
    if(~any(UtilitySU))
        updateFlag = 0; % if UtilitySU has all zeros
        return
    else
        updateFlag = 1;
        if(any(UtilitySU))  % all positive values
            channel = find(UtilitySU, min(UtilitySU));
            B(su, channel) = POperation;
        else    % contains zeros and positive values
            utility = min(UtilitySU(UtilitySU>0)); %minimum positive
            channel = find(UtilitySU, utility);
            B(su, channel) = POperation;
        end
    end
else
    tem2 = find(UtilitySU, min(UtilitySU));
    if (tem2 ~= workingChannel)    %  there is another channel leading to smaller utility
        B(su, tem2) = POperation; % by doing this
        updateFlag = 1;
    else
        updateFlag = 0;
    end    
end








% whiteCat, update channel

% the utlity to switch is: f*p^(0.5)/p, 

% su: index of current secondary user; 

function [B, updateFlag] = update(su, w, B, P, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta)

workingChannel = find(B(su,:)); % which channel is currently used?
n = size(B, 1);
c = size(P, 2);   %   numeracally equal to c, 
Usu = zeros(1, c);
channleBeingUsedByCoLoaction = 0;
realIndex = getRealIndexFromTheExpanded(su, n, w);


for i = 1: c
    %-debug
    if(nnz(B>4) > 0)
       stop =1; 
    end
    %-debug
    % if this channel is not being used by the co-location WBS
    thisChannelIsNotUsedByCOLocationWBSs = channelUsageCoLocationWBSs(i, n, w, su, B);
    if(thisChannelIsNotUsedByCOLocationWBSs)
        B(su, :) = P((realIndex-1)*c + i, :);
        F = (B* B' ~= 0);
        F = F - eye(size(B, 1));
        InterferenceonAll = repmat(Gtilde, w,w).* F * sum(B, 2)*eta + delta/(n/c);
        part1 = InterferenceonAll(su)/(sum(B(su,:))*SUcellRadius^(-pathlossfactor));

        tem = repmat(Gtilde, w, w).*F;
        part2 = sum(tem(:,su)*sum(B(su, :))./(sum(B, 2)*SUcellRadius^(-pathlossfactor))) ; 

        TVInf2Power = (repmat(GtildeAll(n/w+m+1:n/w+m+m, 1:n/w)', w, 1).* (B~=0) * TVpower)./((B+(B==0))*SUcellRadius^(-pathlossfactor));     % B+(B==0) avoid 0 in the demoninator
        part3 = sum(TVInf2Power(:,i));

        
        Usu(i) = part1 + part2 + part3;
        if(Usu(i) == 0)
            stop  =1;
        end
    else
        channleBeingUsedByCoLoaction = i;
    end
end

% delete the element with respect to the channel which is used by the
% co-location WBSs
Usu(Usu == 0 ) = NaN;
 
[MinimumIteration, tem2] = min(Usu);

% if MinimumIteration < PreviousUsu    %  there is another channel leading to smaller utility
%     B(su,:) = P((su-1)*c+tem2, :);
%     updateFlag = 1;
if (tem2 ~= workingChannel)    %  there is another channel leading to smaller utility
    B(su,:) = P((realIndex-1)*c+tem2, :); % by doing this, B(su,:) has only one non-zero element
    updateFlag = 1;
else
    B(su,:) = P((realIndex-1)*c + workingChannel, :);      % still use previous working channel.
    updateFlag = 0;
end    


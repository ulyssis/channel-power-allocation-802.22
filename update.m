% whiteCat, update channel

% the utlity to switch is: f*p^(0.5)/p, 

% su: index of current secondary user; 

function [B, updateFlag] = update(su, B, P, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta)
workingChannel = find(B(su,:)); % which channel is currently used?
n = size(B, 1);
c = size(P, 2);   %   numeracally equal to c, 
% F = (B* B' ~= 0);
% F = F - eye(size(B, 1));
%     InterferenceonAll = Gtilde.* F * sum(B, 2) + delta;
%     part1 = InterferenceonAll(su)/sum(B(su,:));
%     
%     tem = Gtilde.*F;
%     part2 = sum(tem(:,su)*sum(B(su, :))./sum(B, 2));
%     
%     
% 
%         a = sum(B~=0);  % 1 X c
%         TVInf2Power = (GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) * TVpower)./(B+(B==0));     % B+(B==0) avoid 0 in the demoninator
%         
%     part3 = sum(TVInf2Power(:,workingChannel)) / (a(workingChannel)-1);
%    
%     
%     PreviousUsu = part1 + part2 + part3;
%     PreviousChannel = B(su, :); 
%     
%     c = size(P, 2);   %   numeracally equal to c, 
%     Usu = zeros(1 , size(P, 2));    %  store su's utilities which are determined by all different channels
for i = 1: c
        B(su, :) = P((su-1)*c + i, :);
        F = (B* B' ~= 0);
        F = F - eye(size(B, 1));
        InterferenceonAll = Gtilde.* F * sum(B, 2)*eta + delta/(n/c);
        part1 = InterferenceonAll(su)/(sum(B(su,:))*SUcellRadius^(-pathlossfactor));

        tem = Gtilde.*F;
        part2 = sum(tem(:,su)*sum(B(su, :))./(sum(B, 2)*SUcellRadius^(-pathlossfactor))) ; 

        TVInf2Power = (GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) * TVpower)./((B+(B==0))*SUcellRadius^(-pathlossfactor));     % B+(B==0) avoid 0 in the demoninator
        part3 = sum(TVInf2Power(:,i));

        
        Usu(i) = part1 + part2 + part3;
 end

[MinimumIteration, tem2] = min(Usu);

% if MinimumIteration < PreviousUsu    %  there is another channel leading to smaller utility
%     B(su,:) = P((su-1)*c+tem2, :);
%     updateFlag = 1;
if (tem2 ~= workingChannel)    %  there is another channel leading to smaller utility
    B(su,:) = P((su-1)*c+tem2, :); % by doing this, B(su,:) has only one non-zero element
    updateFlag = 1;
else
    B(su,:) = P((su-1)*c + workingChannel, :);      % still use previous working channel.
    updateFlag = 0;
end    


% su: index of current secondary user; 
function [B, updateFlag] = update_self(su, B, P, Gtilde, m, GtildeAll, TVpower, delta)
n = size(B, 1);
F = (B* B' ~= 0);
F = F - eye(size(B, 1));
    InterferenceonAll = Gtilde.* F * sum(B, 2) + sum(GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) * TVpower, 2) + delta;
    part1 = InterferenceonAll(su)/sum(B(su,:));
    
%     tem = Gtilde.*F;
%     part2 = sum(tem(:,su)*sum(B(su, :))./sum(B, 2));
    
%     PreviousUsu = part1 + part2;
    PreviousUsu = part1;
    
    PreviousChannel = B(su, :);
    
    c = size(P, 2);   %   numeracally equal to c, 
Usu = zeros(1 , size(P, 2));    %  store su's utilities which are determined by all different channels
for i = 1: c
        B(su, :) = P((su-1)*c + i, :);
        F = (B* B' ~= 0);
        F = F - eye(size(B, 1));
        InterferenceonAll = Gtilde.* F * sum(B, 2) + sum(GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) *TVpower, 2) + delta;
        part1 = InterferenceonAll(su)/sum(B(su,:));

%         tem = Gtilde.*F;
%         part2 = sum(tem(:,su)*sum(B(su, :))./sum(B, 2));

%         Usu(i) = part1 + part2;
        Usu(i) = part1;
end

[MinimumIteration, tem2] = min(Usu);

if MinimumIteration < PreviousUsu    %  there is another channel leading to smaller utility
    B(su,:) = P((su-1)*c+tem2, :);
    updateFlag = 1;
else
    B(su,:) = PreviousChannel;      % don't change channel.
    updateFlag = 0;
end


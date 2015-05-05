% update channel

% scheme is proved in pimrc_2011 paper
% the utlity to switch is: received power - interference introduced- interference recieved, 

% su: index of current secondary user; 

function [B, updateFlag] = update_distributed_ref(su, B, P, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta, powerLevels)
workingChannel = find(B(su,:)); % which channel is currently used?
workingPower = B(su, workingChannel);
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
Usu = zeros(1, c);
power_onEachChannel=zeros(1, c);
for i = 1: c
%     for j = 1: powerLevels
        B(su, i) = 4 + (P((su-1)*c + i, i)-4)*(j/powerLevels); % min power is 4
        F = (B* B' ~= 0);
        F = F - eye(size(B, 1));
        InterferenceonAll = Gtilde.* F * sum(B, 2)*eta + delta/(n/c);
        part1 = InterferenceonAll(su); % received interference

        tem = Gtilde.*F;
%         part2 = sum(tem(:,su)*sum(B(su, :)));   % caused interference on ohters 
% 
%         part3 = sum(B(su, :))*SUcellRadius^(-pathlossfactor);     % received signal power
        
        if( SUcellRadius^(-pathlossfactor) - sum(tem(:,su)) >0) % the difference between the coefficient of self power and produced interference is positive.
            B(su, i) = P((su-1)*c + i, i);
        else
            B(su, i) = 4;
        end
        
        part2 = sum(tem(:,su)*B(su, i));   % caused interference on ohters 
        part3 = B(su, i)*SUcellRadius^(-pathlossfactor);     % received signal power    
        
        power_onEachChannel(i) =  B(su, i);
        Usu(i) = -part1 - part2 + part3; % minimize this!
        

%     end
end

% find the workingChannel and working Power which result in the minimum Usu
 [tem2, channel]=max(Usu);
%  level=level(1);
%  channel=channel(1);
% update the channel and power uage in matrix B. 
B(su,:) = zeros(1, c);
B(su, channel) = power_onEachChannel(channel);

% if MinimumIteration < PreviousUsu    %  there is another channel leading to smaller utility
%     B(su,:) = P((su-1)*c+tem2, :);
%     updateFlag = 1;
if (channel ~= workingChannel || B(su, channel) ~= workingPower )    %  there is another channel leading to smaller utility
    updateFlag = 1;
else
%     B(su,:) = P((su-1)*c + workingChannel, :);      % still use previous working channel.
    updateFlag = 0;
end    



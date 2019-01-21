% output B through no-regret learning
% sumUtilityNoregret contains at most 1000*16 records
function [B, sumUtilityNoregret, SINRvarianceNoregret] = noregretlearning(seq, B, P, Gtilde, GtildeETsSUs, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, PMiu)
% Reg_i(C x C): regret matrix for each user, update with matrix U
% U(C x C): possible increase of utility, store intermidiate variables
% Pr(C x 1): probability of swithching to new strategy based on regrets
% initialize(Reg(0),A(0)/B(0)) --Pi-> A(1)/B(1) --GlobalNotice--> Upotential(1) --> Reg(1)
            %      |                                                 ^                                                           
            %      |                                                 |
            %      |-------------------------------------------------|
%   (check: B(i,t)=/=B(i,t+1)


n = size(B, 1);
c = size(P, 2);   %   numeracally equal to c, 


% initialize regret matrix and utility matrix for each user

Reg = cell(1, n);
Upotential = cell(1, n);% capatal U
for i = 1: n
    initialReg = rand(c,c);
    Reg{i} = initialReg./(sum(initialReg, 2)*ones(1, c));
    Upotential{i} = zeros(c, c);
end

t = 0;
Bconverge = 0;
        % get the channel index for each user from B:
        I = zeros(1,n);
        for i = 1:n
           [row, col] = find(B(i, :));
           I(i) = col; 
        end
        
    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);
    sumUtilityNoregret = sumUtility;

    
    % record the initial SINR for all SUs
    SINRvarianceNoregret = 0;
    SINR_ETs_previous = [];
    [SINR_ETs_previous, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);        
    
while((~Bconverge) || (t<3))
    Bbackup = B;

   
    for count = 1 : n
        i = seq(count);
       % get utility potential matrix Upotential for i    
        for channel = 1:  c
            Ucurrent = U(i, I(i), B, P, Gtilde, m, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta);
            Upotential{i}(I(i), channel) = Ucurrent - U(i, channel, B, P, Gtilde, m, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta);
        end % notice here, I use (Ucurrent-Uother) to represent the potential, because we are after the minmum U.
     % update each player's Upotential matrix which stores the possible DEcrease of utility by using other strategis.
    
     % update regret matrix
        Reg{i} = (1-1/(t+1))*Reg{i} + Upotential{i};
        miu= 1.2*sum(max(Reg{i}(I(i),:),0));
        if (miu==0)
           miu = 1; 
        end
        
        Pr = zeros(1, c);    % the probability for player i to use a certain strategy
        for strategy=1:c
            if(strategy ~= I(i))
                Pr(strategy) = max(0, Reg{i}(I(i), strategy))/miu;
            else
                Pr(strategy) = 1-sum(max(0, Reg{i}(I(i), :)))/miu;
            end
        end %traverse all strategies of user i
        [row, action] = max(Pr);
        I(i) = action; % choose the trategy with the maximum Pr. if there r several strategies with the same Pr, choose the one with minimum index
        B(i,:) = P((i-1)*c + I(i), :); % renew B
        
        [row, previousAction] = max(Bbackup(i,:));
        if (I(i)~=previousAction)
            [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);
            sumUtilityNoregret(end+1) = sumUtility;
        else
            sumUtilityNoregret(end+1) = sumUtilityNoregret(end);
        end
        
        % calculate the percentage of variance and sum up!
        [SINR_ETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        delta1step = sum(abs((SINR_ETs - SINR_ETs_previous)./SINR_ETs_previous));
        SINRvarianceNoregret = SINRvarianceNoregret + delta1step;        
        SINR_ETs_previous = SINR_ETs;          
    end % traverse all players/users, decide action based on Pr, and form new B.
    
        Bconverge = isequal(B, Bbackup);
                t=t+1;
                if(t==50)
                    Bconverge=1;
                end
end



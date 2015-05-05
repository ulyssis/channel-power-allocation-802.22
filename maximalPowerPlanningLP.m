% Defaultly, primary users utilize channel with the sequence of their index
% #1 PU uses #1 channel.

function [P, linprogWork] = maximalPowerPlanningLP(n, m, infBound, GtildeAll, miniP, maxP)

%             figure(floor(n/100)+1);
%             plot(posSU(1,:), posSU(2,:), 'r+');
%             indexSU = (1:n)';
%             hold on;
%             plot(posPU(1,:), posPU(2,:), 'b<');
%             indexPU = (1:m)';
%             text(posSU(1,:), posSU(2,:), num2str(indexSU));
%             text(posPU(1,:), posPU(2,:), num2str(indexPU));
%             leg3=legend('Secondary users');
%             axis([0 10 0 10]);

c = m;   % #channel == #PU
P = zeros(n*c, c);
pMax = zeros(n, c);    % 

%    f = ones(n, 1)* (-0.01);

lb = miniP * ones(n, 1);
up = maxP * ones(n, 1);
returnValue = zeros(c, 1);

    for channel = 1: m    % as to one channel(or PU), do linear programming to get pMax
        f = -ones(1, n);  % goal is to minmize f'*x, f is a negative coefficient, so minmization becomes maxmazition.
        A = GtildeAll(n+channel, 1:n); % coefficients
        
        % linprog(f,A,b) solves min f'*x such that A*x ≤ b.
        [pMax(:, channel), fval, exitflag] = linprog(f, A, infBound, [], [], lb, up);
        returnValue(channel) = exitflag;
        % pMax --> P
        for i = 1: n
            P((i-1)*c + channel, channel) = pMax(i, channel);
        end
    end
    
    linprogWork = min(returnValue);
    disp(GtildeAll(n+1:n+m, 1:n)*lb);

 % x = linprog(f,A,b,Aeq,beq,lb,ub) defines a set of lower and upper bounds
 % on the design variables, x, so that the solution is always in the range 
 % lb ≤ x ≤ ub. Set Aeq = [] and beq = [] if no equalities exist.   
% to use quadprog to solve

% The optimization problem:
% \sigma_i \sigma_j \sigma_k(q_i^k*\bata_ij^k*h_ji*z) + \sigma_i \sigma_k(q_i^k*x_ik*N_0)
% s.t. x_ij is binary
% \sigma_k x_ik=1;
% p_i^{k,min} <= p_i^k <= p_i^{k,max}
% constrains group 1
% constrains group 2

function [channelUsage, powerLevelFinal] = jointPowerChannelAlloc(n, c, P, Gtilde, delta, minPower)

        condenseP=zeros(n, c);  % P-->condenseP which is n x c
        for i= 1:n
            condenseP(i,:) = sum(P( (i-1)*c+1 : i*c, :), 2);
        end
        
        % number of variables: 3*n*c+2*n^2*c
        % ...x_ik...   ...p_i^k...   \alpha_ij^k...   ...\beta_ij^k...   ...q_i^k        
        %    n*c          n*c           n^2*c             n^2*c         n*c
        
        % matrix H
        % Gtilde
        % decide whether a symmetric matrix: 
        % issym = @(x) isequal(x,x.'): issym(H)
        num_variables = 3*n*c + 2*n^2*c;
        H = zeros(num_variables, num_variables);
        f= zeros(1, num_variables);
        
        for i = 1:n
            for j = 1:n
                for k = 1:c
                    H(2*n*c+n^2*c+((i-1)*n+j-1)*c+k, 2*n*c+2*n^2*c+(i-1)*c+k) = Gtilde(i, j); % q_i^k * \bata_ij^k
                    H(2*n*c+2*n^2*c+(i-1)*c+k, 2*n*c+n^2*c+((i-1)*n+j-1)*c+k) = Gtilde(i, j); % q_i^k * \bata_ij^k, H should be symatric
                    H((i-1)*c+k, 2*n*c+2*n^2*c+(i-1)*c+k) = delta;
                    H(2*n*c+2*n^2*c+(i-1)*c+k, (i-1)*c+k) = delta;
                end
            end
        end
        
        
        % constraints:
        % A*x <= b
        % number of inequations: 6*n^2*c
        % A: 6*n^2*c X 3*n*c+2*n^2*c
        A = [];
        
        b = zeros(6*n^2*c, 1);
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA((j-1)*c+k) =1;
                    oneRowInA((i-1)*c+k) =1;
                    oneRowInA(2*n*c+((i-1)*n+j-1)*c+k) =-1;
                    A = [A; oneRowInA];
                end
            end
        end
        b(1:n^2*c) = ones(n^2*c, 1);
        
        
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA((j-1)*c+k) = -1;
                    oneRowInA((i-1)*c+k) = -1;
                    oneRowInA(2*n*c+((i-1)*n+j-1)*c+k) =2;
                    A = [A; oneRowInA];
                end
            end
        end
        b(n^2*c+1: 2*n^2*c) = zeros(n^2*c, 1);
        
        
        
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA(2*n*c+n^2*c+((i-1)*n+j-1)*c+k) = 1;
                    oneRowInA(2*n*c+((i-1)*n+j-1)*c+k) =-condenseP(j,k);
                    A = [A; oneRowInA];
                end
            end
        end
        b(2*n^2*c+1: 3*n^2*c) = zeros(n^2*c, 1);
        
        
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA(2*n*c+n^2*c+((i-1)*n+j-1)*c+k) = -1;
                    oneRowInA(2*n*c+((i-1)*n+j-1)*c+k) =minPower;
                    A = [A; oneRowInA];
                end
            end
        end
        b(3*n^2*c+1: 4*n^2*c) = zeros(n^2*c, 1);
                
        
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA(2*n*c+n^2*c+((i-1)*n+j-1)*c+k) = 1;
                    oneRowInA(n*c+(j-1)*c+k) = -1;
                    oneRowInA(2*n*c+((i-1)*n+j-1)*c+k) =-minPower;
                    A = [A; oneRowInA];
                end
            end
        end
        b(4*n^2*c+1: 5*n^2*c) = -ones(n^2*c, 1)*minPower;
        
        
        
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA(2*n*c+n^2*c+((i-1)*n+j-1)*c+k) = -1;
                    oneRowInA(n*c+(j-1)*c+k) = 1;
                    oneRowInA(2*n*c+((i-1)*n+j-1)*c+k) =condenseP(j,k);
                    A = [A; oneRowInA];
                end
            end
        end
        
        % the right hand side vulue
        for i = 1:n
            for j= 1:n
                for k =1:c
                    b(5*n^2*c+((i-1)*n+j-1)*c+k) = condenseP(j,k);        
                end
            end
        end

        
        
        %---------
        % Aeq*x = beq
        Aeq = zeros(n,3*n*c+2*n^2*c);
        for i = 1:n
           Aeq(i, (i-1)*c+1: (i-1)*c+c) = ones(1,c); 
        end
        beq = ones(n, 1);
        
        
        %-----
        %lb <= x < = ub
        % as to x, the lb=0; ub =1;
        lb_x = zeros(n*c, 1);
        ub_x = ones(n*c,1);
        
        % as to p, the lb=minPower, ub = xxx
        lb_power = ones(n*c, 1)*minPower;
        ub_power = reshape(condenseP',[1, n*c])';
        
        % \alpha
        lb_alpha = zeros(n^2*c, 1);
        ub_alpha = ones(n^2*c,1);
        
        % \beta
        lb_beta = zeros(n^2*c, 1);
        ub_beta = ones(n^2*c,1);
        
        % q = 1/p
        lb_q = 1./reshape(condenseP',[1, n*c])';
        ub_q = ones(n*c,1)*minPower;        
        
        
        
        lb = vertcat(lb_x, lb_power, lb_alpha, lb_beta, lb_q);
        ub = vertcat(ub_x, ub_power, ub_alpha, ub_beta, ub_q);
        
        
        x = quadprog(H,f,A,b,Aeq,beq,lb,ub);
        channelUsage = x(1: n*c);
        powerLevelFinal = x(n*c+1: 2*n*c);

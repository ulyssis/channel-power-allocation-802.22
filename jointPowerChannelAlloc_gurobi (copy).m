function [channelUsage, powerLevelFinal] = jointPowerChannelAlloc_gurobi(n, c, P, Gtilde, delta, minPower)



% The optimization problem:
% min    \sigma_i \sigma_j \sigma_k(q_i^k*\bata_ij^k*h_ji*z) + \sigma_i \sigma_k(q_i^k*x_ik*N_0)
% s.t. x_ij is binary
% \sigma_k x_ik=1;
% p_i^{k,min} <= p_i^k <= p_i^{k,max}
% constrains group 1
% constrains group 2

clear model;
clear params;



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
                    H(2*n*c + n^2*c + ((i-1)*n+j-1)*c + k, 2*n*c + 2*n^2*c + (i-1)*c + k) = Gtilde(i, j); %  \bata_ij^k * q_i^k
                    H(2*n*c + 2*n^2*c + (i-1)*c + k, 2*n*c + n^2*c + ((i-1)*n+j-1)*c + k) = Gtilde(i, j); % q_i^k * \bata_ij^k, H should be symatric
                    H(2*n*c+2*n^2*c+(i-1)*c+k, (i-1)*c+k) = delta;
                    H((i-1)*c+k, 2*n*c+2*n^2*c+(i-1)*c+k) = delta;
                end
            end
        end

model.Q = sparse(H);        % quadratic objective matrix, must be sparse
model.obj = f;      % .obj: the linear objective vetor in the optimized function in the problem. must be dense.



        % constraints:
        % A*x <= b
        % number of inequations: 6*n^2*c
        % A: 6*n^2*c X 3*n*c+2*n^2*c
        A = [];
        
        % rhs
        %b = zeros(6*n^2*c, 1);
        b=[];

%% Linearization method 1
%         % x_jk + x_ik -\alpha_ij^k <=1
%         for i = 1:n
%             for j= 1:n
%                 for k = 1:c
%                     oneRowInA = zeros(1, 3*n*c+2*n^2*c);
%                     oneRowInA((j-1)*c+k) =1;
%                     oneRowInA((i-1)*c+k) =1;
%                     oneRowInA(2*n*c + hash_script2index(i, j, k, n, c)) =-1;
%                     A = [A; oneRowInA];
%                 end
%             end
%         end
%         bineq =  ones(n^2*c, 1);
%         %b(1:n^2*c) = ones(n^2*c, 1);
%         b=[b; bineq];
%         
%         % -x_jk - x_ik +\alpha_ij^k <=0
%         for i = 1:n
%             for j= 1:n
%                 for k = 1:c
%                     oneRowInA = zeros(1, 3*n*c+2*n^2*c);
%                     oneRowInA((j-1)*c+k) = -1;
%                     oneRowInA((i-1)*c+k) = -1;
%                     oneRowInA(2*n*c+ hash_script2index(i, j, k, n, c) ) =2;
%                     A = [A; oneRowInA];
%                 end
%             end
%         end
%         %b(n^2*c+1: 2*n^2*c) = zeros(n^2*c, 1);
%         bineq =  zeros(n^2*c, 1);
%         b=[b; bineq];
%% Linearization method 1
        
%% Linearization method 2
        % -x_ik + \alpha_ij^k <0
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA((i-1)*c+k) =-1;
                    oneRowInA(2*n*c + hash_script2index(i, j, k, n, c)) =1;
                    A = [A; oneRowInA];
                end
            end
        end
        bineq =  zeros(n^2*c, 1);
        %b(1:n^2*c) = ones(n^2*c, 1);
        b=[b; bineq];

        % -x_jk + \alpha_ij^k <0
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA((j-1)*c+k) =-1;
                    oneRowInA(2*n*c + hash_script2index(i, j, k, n, c)) =1;
                    A = [A; oneRowInA];
                end
            end
        end
        bineq =  zeros(n^2*c, 1);
        %b(1:n^2*c) = ones(n^2*c, 1);
        b=[b; bineq];

        % x_ik + x_jk -1 - \alpha_ij^k <=0
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA((i-1)*c+k) =1;
                    oneRowInA((j-1)*c+k) =1;
                    oneRowInA(2*n*c + hash_script2index(i, j, k, n, c)) =-1;
                    A = [A; oneRowInA];
                end
            end
        end
        bineq =  ones(n^2*c, 1);
        %b(1:n^2*c) = ones(n^2*c, 1);
        b=[b; bineq];
%% end of Linearization method 2

        
        
        %\beta_ij^k - p_j^{k,max}*\alpha_{ij}^k <=0
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA(2*n*c + n^2*c + hash_script2index(i, j, k, n, c)) = 1;
                    oneRowInA(2*n*c + hash_script2index(i, j, k, n, c)) =-condenseP(j,k);
                    A = [A; oneRowInA];
                end
            end
        end
%         b(2*n^2*c+1: 3*n^2*c) = zeros(n^2*c, 1);
        bineq =  zeros(n^2*c, 1);
        b=[b; bineq];
        
        % p_j^{k,min} * \alpha_{ij}^k - \beta_ij^k  <=0
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA(2*n*c + n^2*c + hash_script2index(i, j, k, n, c)) = -1;
                    oneRowInA(2*n*c + hash_script2index(i, j, k, n, c)) =minPower;
                    A = [A; oneRowInA];
                end
            end
        end
        %b(3*n^2*c+1: 4*n^2*c) = zeros(n^2*c, 1);
        bineq =  zeros(n^2*c, 1);
        b=[b; bineq];
                
        % \beta_ij^k - p_j^k - p_j^{k,min} * \alpha_{ij}^k <= -p_j^{k,min}
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA(2*n*c + n^2*c + hash_script2index(i, j, k, n, c)) = 1;
                    oneRowInA(n*c + (j-1)*c+k) = -1;
                    oneRowInA(2*n*c + hash_script2index(i, j, k, n, c)) =-minPower;
                    A = [A; oneRowInA];
                end
            end
        end
%         b(4*n^2*c+1: 5*n^2*c) = -ones(n^2*c, 1)*minPower;
        bineq =  -ones(n^2*c, 1)*minPower;
        b=[b; bineq];
        
        
        % -\beta_ij^k + p_j^k + p_j^{k,max} * \alpha_{ij}^k <= p_j^{k,max}        
        for i = 1:n
            for j= 1:n
                for k = 1:c
                    oneRowInA = zeros(1, 3*n*c+2*n^2*c);
                    oneRowInA(2*n*c + n^2*c + hash_script2index(i, j, k, n, c)) = -1;
                    oneRowInA(n*c + (j-1)*c+k) = 1;
                    oneRowInA(2*n*c + hash_script2index(i, j, k, n, c)) =condenseP(j,k);
                    A = [A; oneRowInA];
                end
            end
        end
        
        b_length = size(b, 1);
        % a chunk of right hand side vulues
        for i = 1:n
            for j= 1:n
                for k =1:c
                    b(b_length + hash_script2index(i, j, k, n, c)) = condenseP(j,k);        
                end
            end
        end

        
        % linear constraints: equality
        % Aeq*x = beq
        % number of constraints is n
        Aeq = zeros(n,3*n*c+2*n^2*c);
        for i = 1:n
           Aeq(i, (i-1)*c+1: (i-1)*c+c) = ones(1,c);
        end
        A = [A; Aeq];
        beq = ones(n, 1);        
        b =[b;beq];
        
model.A = sparse(A);
model.rhs = b;
        
charSense = [];
for i = 1: 7*n^2*c
   charSense =[charSense, '<'];  
end
for i = 1: n
   charSense =[charSense, '=']; 
end
model.sense = charSense;
        
%%
% %% add quadratic contraint
% %       p_i^k * q_i^k = 1
% Qc = zeros(num_variables, num_variables);
% for i = 1:n
%     for k = 1:c
%         Qc(n*c + (i-1)*c + k, 2*n*c+2*n^2*c+(i-1)*c + k) = 1;
%         Qc(2*n*c+2*n^2*c+(i-1)*c + k, n*c + (i-1)*c + k) = 1;
%     end
% end
% 
% model.quadcon.Qc = sparse(Qc);
% 
% Qcrhs = ones(n*c, 1);
% model.quadcon.rhs = Qcrhs;
% model.quadcon.q = f;
% 
% Qcsense = [];
% for i=1:n*c
%     Qcsense =[Qcsense, '='];
% end
% model.quadcon.sense = Qcsense;
%%


        
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
        ub_beta = 40*ones(n^2*c,1);
        
        % q = 1/p
        lb_q = 1./reshape(condenseP',[1, n*c])';
        ub_q = 1./(ones(n*c,1)*minPower);        
        
        
        
        lb = vertcat(lb_x, lb_power, lb_alpha, lb_beta, lb_q);
        ub = vertcat(ub_x, ub_power, ub_alpha, ub_beta, ub_q);
        
model.lb = lb';        
model.ub = ub';        

% ’C’ (continuous), ’B’ (binary), ’I’ (integer), ’S’ (semi-continuous), or ’N’ (semi-integer)
vtype =[];
for i = 1: n*c
   vtype =[vtype, 'B']; 
end

for i = 1: n*c
   vtype =[vtype, 'C']; 
end

for i = 1: n^2*c
   vtype =[vtype, 'B']; 
end

for i = 1: n^2*c+n*c
   vtype =[vtype, 'C']; 
end

model.vtype = vtype;



%%model.start =nan;
start_vector_x_ik = zeros(1, n*c);
% % random initialization
% for i =1:n
%     offset =floor(rand(1)*c) + 1;
%     start_vector_x_ik((i-1)*c + offset) =1;
% end
% fixed initialization    
start_vector_x_ik(1:c:length(start_vector_x_ik)) =1;


start_vector_p = ub_power';

start_vector_alpha = [];
for i=1:n
    for j = 1:n
        for k = 1:c
        start_vector_alpha = [start_vector_alpha, start_vector_x_ik((i-1)*c+k)*start_vector_x_ik((j-1)*c+k)];
        end
    end
end

start_vector_beta = [];
for i=1:n
    for j = 1:n
        for k = 1:c
        start_vector_beta = [start_vector_beta, start_vector_alpha(hash_script2index(i, j, k, n, c))*start_vector_p((j-1)*c+k)];
        end
    end
end

start_vector_q = 1./start_vector_p;
start = [start_vector_x_ik, start_vector_p, start_vector_alpha, start_vector_beta, start_vector_q];



model.start =start;



%         x = quadprog(H,f,A,b,Aeq,beq,lb,ub);
%         channelUsage = x(1: n*c);
%         powerLevelFinal = x(n*c+1: 2*n*c);

model.modelsense = 'min';

gurobi_write(model, 'jointPC.mps');


% model = gurobi_read('jointPC.mps');
% clear params;
% params.
fprintf('there are %f variables in the joint power channel optimization.\n\n\n', 3*n*c + 2*n^2*c );

% params.method = 1;
% results = gurobi(model, params);
results = gurobi(model);

disp(results);
fprintf('\n.status %e \n .runtime %e \n', results.status, results.runtime);

fprintf('x_ij:\n');
tmp_list = results.x(1:n*c);
% for j=1:2*n*c
% fprintf('%e\n', results.x(j));
% end
reshape(tmp_list, c, n)'  % one row corresponds one node

fprintf('p_i^k:\n');
tmp_list = results.x(n*c+1: 2*n*c);
reshape(tmp_list, c, n)

fprintf('\alpha_ij^k:\n');
tmp_list1 = results.x(2*n*c+1: 2*n*c+n^2*c);
for k=1:c
    tmp_list = tmp_list1(k: c: n^2*c);
    reshape(tmp_list, n, n)
    fprintf('-----\n');
end



fprintf('q_i^k:\n');
tmp_list = results.x(2*n*c+2*n^2*c+1: 3*n*c+2*n^2*c);
reshape(tmp_list, n, c)


fprintf('\n------\n');

fprintf('Obj: %e\n', results.objval);


function [B] = GUROBI_ECC(n, c, w, P, Gtilde, delta)      

        condenseP=zeros(n, c);  % P-->condenseP which is n x c: p_i^k, i\in N, k\in C
        for i= 1:n
            condenseP(i,:) = sum(P( (i-1)*c+1 : i*c, :), 2);
        end


       % X:
       % x1 x2 x3 ... xn    x1 x2 x3 ... xn     x1 x2 x3 ... xn
       % ----channel 1--    ---channel 2---     ----channel 3-- 
        % model.Q
        % h_{i,j}*z_{i,j}*p_{j,k}/p_{i,k}
        % h_{i,j}*z_{i,j} is element of Gtilde.

        arrayH = zeros(n, n);
        arrayH = repmat (arrayH, c);
        for i = 1: c
            newpart_nominator = ones(n, 1) * condenseP(:, i)';
            newpart_dominator = condenseP(:, i) * ones(1, n);
            newpart = newpart_nominator./newpart_dominator;
            newpart = newpart.*Gtilde;
            arrayH(1+(i-1)*n : i*n, 1+(i-1)*n : i*n) = newpart;
        end
        arrayH = arrayH - arrayH.*eye(n*c);

        
        % model.obj
        % N_0/p_{i,k}
        %NoisePowerRatio = delta./(condenseP*SUcellRadius^(-pathlossfactor));
        NoisePowerRatio = delta./(condenseP);
        NoisePowerRatioInOneRow =[];
        for i=1:c
            NoisePowerRatioInOneRow = [NoisePowerRatioInOneRow, NoisePowerRatio(:, i)'];
        end
        
        % model.A, parameters in constraints
        % n x (n*c)
        % assuming n = 4, c = 2, then model.A is,
        % 1 0 0 0 1 0 0 0
        % 0 1 0 0 0 1 0 0
        % 0 0 1 0 0 0 1 0
        % 0 0 0 1 0 0 0 1
        A = [];
        for i= 1:n
            row = zeros(1, n*c);
            for j = 1:c
                row(i + (j-1)*n ) = 1;
            end

            A = [A; row];
        end

        arrayH = arrayH.*(10e+10);
        NoisePowerRatioInOneRow = NoisePowerRatioInOneRow.*(10e+10);

        optimizaionModel1.Q = sparse(arrayH);
        optimizaionModel1.obj = NoisePowerRatioInOneRow;
        
        optimizaionModel1.A = sparse(A);
        optimizaionModel1.rhs = ones(1, n)*w;
        optimizaionModel1.sense = '=';
        optimizaionModel1.vtype = 'B';
        %optimizaionModel1.modelsense = 'min';
        opt1Results = gurobi(optimizaionModel1);
        
        %assignin('base', 'results', opt1Results);
        resultX = zeros(n, c);
        reversedResultX = zeros(c, n);
        for i=1:c
            reversedResultX(i, :) = opt1Results.x((i-1)*n+1 : i*n);
        end
        resultX = reversedResultX';
        B= resultX.*condenseP;
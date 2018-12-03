function[opt1Results] = solveSchmeIIwithGORUBI(Gtilde, n, c, GtildeAll, delta, PMiu, POperation, infBound)
                                               
delta = delta * (1e+10);
Gtilde = Gtilde * (1e+10);
GtildeAll = GtildeAll * (1e+10);
infBound = infBound * (1e+10);

   % variables, the number of variables is n*c + 2n^2c +2n
   % _________________________________________________________________
   % variable | x_i^k        \alpha_i^k       \beta_i^k    y_i    z_i
   % number   | n*c          n^2*c            n^2*c         n      n
   % -----------------------------------------------------------------
       
   %% model.Q, the parameters of qradratics in the objective function
   arrayForQ = zeros(n*c + 2*n^2*c +2*n, n*c + 2*n^2*c +2*n);
       for i = 1:n
        for j = i:n
            if(j~=i)
                for k = 1:c
                    rowIndex1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;  % beta_{ij}^k
                    rowIndex2 = (k-1)*n + i;                            % x_i^k
                    columnIndex1 = n*c + 2*n^2*c + n + i;               % z_i
                    arrayForQ (rowIndex1, columnIndex1) = GtildeAll(j, i)*(PMiu-POperation);
                    arrayForQ (rowIndex2, columnIndex1) = delta*(PMiu-POperation);
                end
            end
        end
       end
       
       
       %% model.obj, the parameters of linear parts in the objective function
       vectorForObj = zeros(1, n*c + 2*n^2*c +2*n);
       
       for i = 1:n
        for j = i:n
            if(j~=i)
                for k = 1:c
                    index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j; % beta
                    index2 = (k-1)*n + i;                           % x_i^k
                    vectorForObj (index1) = GtildeAll(j, i) * POperation;
                    vectorForObj (index2) = delta * POperation;
                end
            end
        end
       end
       
       
       %% model.A, the parameters for the linear constraints
       % contraits: there are 6*n^2*c + 5*n linear inequalities/equalities in the constraitns
       % and n*c + 2*n^2*c + 2*n variables.
       arrayForA = zeros(6*n^2*c + 5*n, n*c + 2*n^2*c +2*n);
       % model.rhs
       modelRHS = [];
       modelSENSE = [];
       
       ConstraintGroupStartsIndex = 1;
       indexOfLinearConstraints = ConstraintGroupStartsIndex;
            % 1st linear inequality group
               for i = 1:n
                for j = i+1:n
                        for k = 1:c
                            % index of the linear inequality
                            indexInequality = indexOfLinearConstraints;
                            index1 = (k-1)*n + i;
                            index2 = (k-1)*n + j;
                            index3 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = 1;
                            arrayForA(indexInequality, index2) = 1;
                            arrayForA(indexInequality, index3) = -1;
                            indexOfLinearConstraints = indexOfLinearConstraints + 1;                            
                        end
                end
               end
               
               rhs1 = ones(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               modelRHS = [modelRHS rhs1];
               
               senseVector1 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
                senseVector1(:) = {'<'};
                sense1 = char(senseVector1);
                modelSENSE = [modelSENSE sense1'];
        
        
             % 2nd linear inequality group
                ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                for j = i+1:n
                        for k = 1:c
                            % index of the linear inequality
                            indexInequality = indexOfLinearConstraints;
                            index1 = (k-1)*n + i;
                            index2 = (k-1)*n + j;
                            index3 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = -1;
                            arrayForA(indexInequality, index2) = -1;
                            arrayForA(indexInequality, index3) = 2;
                            indexOfLinearConstraints = indexOfLinearConstraints + 1;                            
                        end
                end
               end
               rhs2 = zeros(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               modelRHS = [modelRHS rhs2];
               
               senseVector2 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector2(:) = {'<'};
               sense2 = char(senseVector2);
               modelSENSE = [modelSENSE sense2'];               
               
             % 3rd linear inequality group
                ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                for j = i+1:n
                        for k = 1:c
                            % index of the linear inequality
                            indexInequality = indexOfLinearConstraints;
                            index2 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = 1;
                            arrayForA(indexInequality, index2) = -POperation;
                            indexOfLinearConstraints = indexOfLinearConstraints + 1;                            
                        end
                end
               end
               rhs3 = zeros(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               modelRHS = [modelRHS rhs3];
               
               senseVector3 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector3(:) = {'<'};
               sense3 = char(senseVector3);
               modelSENSE = [modelSENSE sense3'];  
               
             % 4th linear inequality group
             ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                for j = i+1:n
                        for k = 1:c
                            % index of the linear inequality
                            indexInequality = indexOfLinearConstraints;
                            index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
                            index2 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = -1;
                            arrayForA(indexInequality, index2) = PMiu;
                            indexOfLinearConstraints = indexOfLinearConstraints + 1;                            
                        end
                end
               end       
               rhs4 = zeros(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
              modelRHS = [modelRHS rhs4];
              senseVector4 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector4(:) = {'<'};
               sense4 = char(senseVector4);
               modelSENSE = [modelSENSE sense4'];  
              
              
             % 5th linear inequality group
             ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                for j = i+1:n
                        for k = 1:c
                            % index of the linear inequality
                            indexInequality = indexOfLinearConstraints;
                            index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
                            index2 = n*c + 2*n^2*c + j;
                            index3 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = 1;
                            arrayForA(indexInequality, index2) = POperation - PMiu;
                            arrayForA(indexInequality, index3) = - PMiu;
                            indexOfLinearConstraints = indexOfLinearConstraints + 1;                            
                        end
                end
               end
               rhs5 = ones(1, indexOfLinearConstraints - ConstraintGroupStartsIndex)*(POperation - PMiu);
               modelRHS = [modelRHS rhs5];
               senseVector5 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector5(:) = {'<'};
               sense5 = char(senseVector5);
               modelSENSE = [modelSENSE sense5'];                
               
               
%              % 6th linear inequality group
%              ConstraintGroupStartsIndex = indexOfLinearConstraints;
%                for i = 1:n
%                 for j = i+1:n
%                         for k = 1:c
%                             % index of the linear inequality
%                             indexInequality = indexOfLinearConstraints;
%                             index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
%                             index2 = n*c + 2*n^2*c + j;
%                             index3 = n*c + (k-1)*n^2 + (i-1)*n + j;
%                             arrayForA(indexInequality, index1) = -1;
%                             arrayForA(indexInequality, index2) = -POperation + PMiu;
%                             arrayForA(indexInequality, index3) = POperation;
%                             indexOfLinearConstraints = indexOfLinearConstraints + 1;                            
%                         end
%                 end
%                end    
%                rhs6 = zeros(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
%                modelRHS = [modelRHS rhs6];
%                
%                senseVector6 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
%                senseVector6(:) = {'<'};
%                sense6 = char(senseVector6);
%                modelSENSE = [modelSENSE sense6'];                
               
               
               
             % quadratic equality relaxation p_to_binary_rex_1
             % 7th linear inequality group
             ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                % index of the linear inequalitie
                indexInequality = indexOfLinearConstraints;
                index1 = n*c + 2*n^2*c + i;
                index2 = n*c + 2*n^2*c + n + i;
                arrayForA(indexInequality, index1) = PMiu;
                arrayForA(indexInequality, index2) = POperation;
                indexOfLinearConstraints = indexOfLinearConstraints + 1;
               end 
               rhs7 = ones(1, indexOfLinearConstraints - ConstraintGroupStartsIndex) * (POperation^2-1)/(POperation - PMiu);
                modelRHS = [modelRHS rhs7];

               senseVector7 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector7(:) = {'<'};
               sense7 = char(senseVector7);
               modelSENSE = [modelSENSE sense7'];   
               
               
             %  quadratic equality relaxation, p_to_binary_rex_2
             % 8th linear inequality group
               ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                % index of the linear inequalitie
                indexInequality = indexOfLinearConstraints;
                index1 = n*c + 2*n^2*c + i;
                index2 = n*c + 2*n^2*c + n + i;
                arrayForA(indexInequality, index1) = POperation;
                arrayForA(indexInequality, index2) = PMiu;
                indexOfLinearConstraints = indexOfLinearConstraints + 1;
               end      
               rhs8 = ones(1, indexOfLinearConstraints - ConstraintGroupStartsIndex) * (POperation^2-1)/(POperation - PMiu);
                modelRHS = [modelRHS rhs8];
                
               senseVector8 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector8(:) = {'<'};
               sense8 = char(senseVector8);
               modelSENSE = [modelSENSE sense8'];           

               
               % quadratic equality relaxation, p_to_binary_rex_3
               % 9th linear inequality group
               ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                % index of the linear inequalitie
                indexInequality = indexOfLinearConstraints;
                index1 = n*c + 2*n^2*c + i;
                index2 = n*c + 2*n^2*c + n + i;
                arrayForA(indexInequality, index1) = -1;
                arrayForA(indexInequality, index2) = -1;
                indexOfLinearConstraints = indexOfLinearConstraints + 1;
               end      
               rhs9 = ones(1, indexOfLinearConstraints - ConstraintGroupStartsIndex) * (1+PMiu^2-2*PMiu*POperation)/(PMiu*POperation - PMiu^2);
               modelRHS = [modelRHS rhs9];
                
               senseVector9 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector9(:) = {'<'};
               sense9 = char(senseVector9);
               modelSENSE = [modelSENSE sense9']; 
               
               
               % quadratic equality relaxation, p_to_binary_rex_4
               % 10th linear inequality group
               ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                % index of the linear inequalitie
                indexInequality = indexOfLinearConstraints;
                index1 = n*c + 2*n^2*c + i;
                index2 = n*c + 2*n^2*c + n + i;
                arrayForA(indexInequality, index1) = -1;
                arrayForA(indexInequality, index2) = -1;
                indexOfLinearConstraints = indexOfLinearConstraints + 1;
               end      
               rhs10 = ones(1, indexOfLinearConstraints - ConstraintGroupStartsIndex) * (1- POperation^2)/(POperation^2 - PMiu*POperation);
               modelRHS = [modelRHS rhs10];

               senseVector10 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector10(:) = {'<'};
               sense10 = char(senseVector10);
               modelSENSE = [modelSENSE sense10'];      
               
               
              % 11th linear Equality group
               ConstraintGroupStartsIndex = indexOfLinearConstraints;
               for i = 1:n
                indexInequality = indexOfLinearConstraints;
                for k = 1:c
                    arrayForA(indexInequality, (k-1)*n + i) = 1;
                end
                indexOfLinearConstraints = indexOfLinearConstraints + 1;

               end                
               rhs11 = ones(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
                modelRHS = [modelRHS rhs11];
                
               senseVector11 = cell(1, indexOfLinearConstraints - ConstraintGroupStartsIndex);
               senseVector11(:) = {'='};
               sense11 = char(senseVector11);
               modelSENSE = [modelSENSE sense11'];                    
                
         arrayForA_2 = arrayForA(1:indexInequality, :);
               
        %% solve model
        model.Q   = sparse(arrayForQ);

        model.obj = vectorForObj;

        model.A = sparse(arrayForA_2);
        
        model.rhs = modelRHS;
        
        model.sense = modelSENSE;

        model.vtype = 'B';
        
        %% solve model.quadcon
        % The quadratic constraint
        
        for k = 1:c
            matrix = zeros(n*c + 2*n^2*c + 2*n, n*c + 2*n^2*c + 2*n);
            vector = zeros(n*c + 2*n^2*c + 2*n, 1);
            for i= 1:n
                    matrix((k-1)*n+i, n*c + 2*n^2*c + i) = (PMiu - POperation)*GtildeAll(n+k, i);
                    %vector((k-1)*n +i : k*n) = POperation * GtildeAll(i, n+k);
                    vector((k-1)*n +i) = POperation * GtildeAll(n+k, i);
            end
            model.quadcon(k).Qc = sparse(matrix);
            model.quadcon(k).q = vector;
            model.quadcon(k).rhs = infBound;
            model.quadcon(k).name = 'rot_cone';
        end

                

       opt1Results = gurobi(model);
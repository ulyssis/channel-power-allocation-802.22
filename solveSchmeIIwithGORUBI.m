function[opt1Results] = solveSchmeIIwithGORUBI(Gtilde, n, c, GtildeAll, delta, PMiu, POperation, infBound)
                                               
            
   % variables, the number of variables is n*c + 2n^2c +2n
   % _________________________________________________________________
   % variable | x_i^k        \alpha_i^k       \beta_i^k    y_i    z_i
   % number   | n*c          n^2*c            n^2*c         n      n
   % -----------------------------------------------------------------
       
   %% model.Q
   arrayForQ = zeros(n*c + 2*n^2*c +2*n, n*c + 2*n^2*c +2*n);
       for i = 1:n
        for j = i:n
            if(j~=i)
                for k = i:c
                    rowIndex1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
                    columnIndex1 = n*c + 2*n^2*c + n + i; 
                    rowIndex2 = (i-1)*n + i;
                    arrayForQ (rowIndex1, columnIndex1) = Gtilde(j, i)*(PMiu-POperation);
                    arrayForQ (rowIndex2, columnIndex1) = delta*(PMiu-POperation);
                end
            end
        end
       end
       
       
       %% model.obj
       vectorForObj = zeros(1, n*c + 2*n^2*c +2*n);
       
       for i = 1:n
        for j = i:n
            if(j~=i)
                for k = i:c
                    index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
                    index2 = (i-1)*n + i;
                    vectorForObj (index1) = Gtilde(j, i) * POperation;
                    vectorForObj (index2) = delta * POperation;
                end
            end
        end
       end
        
       %% model.A
       % contraits: there are 6*n^2*c+3*n linear inequalities in the constraitns
       arrayForA = zeros(6*n^2*c + 3*n, n*c + 2*n^2*c +2*n);
            % 1st linear inequality group
               for i = 1:n
                for j = i:n
                    if(j~=i)
                        for k = i:c
                            % index of the linear inequalitie
                            indexInequality = (i-1)*n + (j-1)*n + k;
                            index1 = (k-1)*n + i;
                            index2 = (k-1)*n + j;
                            index3 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = 1;
                            arrayForA(indexInequality, index2) = 1;
                            arrayForA(indexInequality, index3) = -1;
                        end
                    end
                end
               end
               
               rhs1 = ones(1, n^2*c);
               
             % 2nd linear inequality group
               for i = 1:n
                for j = i:n
                    if(j~=i)
                        for k = i:c
                            % index of the linear inequalitie
                            indexInequality = n^2*c + (i-1)*n + (j-1)*n + k;
                            index1 = (k-1)*n + i;
                            index2 = (k-1)*n + j;
                            index3 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = -1;
                            arrayForA(indexInequality, index2) = -1;
                            arrayForA(indexInequality, index3) = 2;
                        end
                    end
                end
               end
               rhs2 = zeros(1, n^2*c);
               
             % 3rd linear inequality group
               for i = 1:n
                for j = i:n
                    if(j~=i)
                        for k = i:c
                            % index of the linear inequalitie
                            indexInequality = 2*n^2*c + (i-1)*n + (j-1)*n + k;
                            index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
                            index2 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = 1;
                            arrayForA(indexInequality, index2) = -POperation;
                        end
                    end
                end
               end
               rhs3 = zeros(1, n^2*c);
               
             % 4th linear inequality group
               for i = 1:n
                for j = i:n
                    if(j~=i)
                        for k = i:c
                            % index of the linear inequalitie
                            indexInequality = 3*n^2*c + (i-1)*n + (j-1)*n + k;
                            index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
                            index2 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = -1;
                            arrayForA(indexInequality, index2) = PMiu;
                        end
                    end
                end
               end       
               rhs4 = zeros(1, n^2*c);
              
             % 5th linear inequality group
               for i = 1:n
                for j = i:n
                    if(j~=i)
                        for k = i:c
                            % index of the linear inequalitie
                            indexInequality = 4*n^2*c + (i-1)*n + (j-1)*n + k;
                            index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
                            index2 = n*c + 2*n^2*c + j;
                            index3 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = 1;
                            arrayForA(indexInequality, index2) = POperation - PMiu;
                            arrayForA(indexInequality, index3) = - PMiu;
                        end
                    end
                end
               end
               rhs5 = ones(1, n^2*c)*(POperation - PMiu);
               
             % 6th linear inequality group
               for i = 1:n
                for j = i:n
                    if(j~=i)
                        for k = i:c
                            % index of the linear inequalitie
                            indexInequality = 5*n^2*c + (i-1)*n + (j-1)*n + k;
                            index1 = n*c + n^2*c + (k-1)*n^2 + (i-1)*n + j;
                            index2 = n*c + 2*n^2*c + j;
                            index3 = n*c + (k-1)*n^2 + (i-1)*n + j;
                            arrayForA(indexInequality, index1) = -1;
                            arrayForA(indexInequality, index2) = -POperation + PMiu;
                            arrayForA(indexInequality, index3) = POperation;
                        end
                    end
                end
               end    
               rhs6 = zeros(1, n^2*c);
               
               
             % 7th linear inequality group
               for i = 1:n
                % index of the linear inequalitie
                indexInequality = 6*n^2*c + i;
                index1 = n*c + 2*n^2*c + i;
                index2 = n*c + 2*n^2*c + n + i;
                arrayForA(indexInequality, index1) = POperation;
                arrayForA(indexInequality, index2) = PMiu;                
               end 
               rhs7 = ones(1, n) * (POperation^2-1)/(POperation - PMiu);
               
               
             % 8th linear inequality group
               for i = 1:n
                % index of the linear inequalitie
                indexInequality = 6*n^2*c + n + i;
                index1 = n*c + 2*n^2*c + i;
                index2 = n*c + 2*n^2*c + n + i;
                arrayForA(indexInequality, index1) = PMiu;
                arrayForA(indexInequality, index2) = POperation;                
               end      
               rhs8 = ones(1, n) * (POperation^2-1)/(POperation - PMiu);

               
              % 9th linear equality group
               for i = 1:n
                % index of the linear inequalitie
                indexInequality = 6*n^2*c + 2*n + i;
                for k = 1:c
                    arrayForA(indexInequality, (k-1)*n + i) = 1;
                end
               end                
               rhs9 = ones(1, n);
               
         
        %% solve model
        model.Q   = sparse(arrayForQ);

        model.obj = vectorForObj;

        model.A = sparse(arrayForA);

        model.rhs = [rhs1 rhs2 rhs3 rhs4 rhs5 rhs6 rhs7 rhs8 rhs9];

        senseVector1 = cell(1, (6*n^2*c + 2*n));
        senseVector1(:) = {'<'};
        sense1 = char(senseVector1);
        senseVector2 = cell(1, n);
        senseVector2(:) = {'='};
        sense2 = char(senseVector2);

        model.sense = [sense1' sense2'];

        model.vtype = 'B';
        
        %% solve model.quadcon
        
        % The quadratic in the constraint
        
        for k = 1:c
            matrix = zeros(n*c + 2*n^2*c + 2*n, n*c + 2*n^2*c + 2*n);
            vector = zeros(n*c + 2*n^2*c + 2*n, 1);
            for i= 1:n
                    matrix((k-1)*n+i, n*c + 2*n^2*c + i) = (PMiu - POperation)*GtildeAll(i, n+k);
                    vector((k-1)*n +i : k*n) = POperation * GtildeAll(i, n+k);
            end
            model.quadcon(k).Qc = sparse(matrix);
            model.quadcon(k).q = vector;
            model.quadcon(k).rhs = infBound;
            model.quadcon(k).name = 'rot_cone';
        end

                
    
       opt1Results = gurobi(model);
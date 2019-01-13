       
function [utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, B_random, B_cat, B_case, B_optimization, B_noregret, B_PotentialGame, RetGUROBI_FCC] = ...
    runSchemes(run, P, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, ...
    PMiu, POperation, infBound, RetGUROBI_FCC)

% function [RetGUROBI_FCC] = runSchemes(run, P, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, ...
%     schemeIIEnabled, PMiu, POperation, infBound, RetGUROBI_FCC)

seq = randperm(n);
SchemeIISolutionFlag = 0;
%% scheme II
if(schemeIIEnabled)
    %% centralized
    opt2Results = solveSchmeIIwithGORUBI(Gtilde, n, c, GtildeAll, delta, PMiu, POperation, infBound);

        % -- to record how many times GUROBI returns solution.
        resultStatus = strcmp(opt2Results.status, 'OPTIMAL');
        if(resultStatus)
            RetGUROBI_FCC(run) = 1;
            SchemeIISolutionFlag = 1;
                        
            %% record how many times GUROBI returns solution
            ResultXSchemeII = zeros(n, c);
            for k=1:c
                ResultXSchemeII(:, k) = opt2Results.x((k-1)*n + 1 : k*n)';
            end
            %% get B_scheme2Centralized
            % the n auxiliary variables
            ResultYSchemeII = opt2Results.x(n*c + 2*n^2*c +1: n*c + 2*n^2*c +n);   
            if(nnz(ResultYSchemeII) > 0)
                stop =1;
            end
            % the n auxiliary variables
            ResultZSchemeII = opt2Results.x(n*c + 2*n^2*c +n+1: n*c + 2*n^2*c + 2*n);

            synthesisSchemeII = zeros(1, n);

            for i = 1:n
                if(ResultYSchemeII(i) ~= ResultZSchemeII(i))
                    synthesisSchemeII(i) = ResultYSchemeII(i)*PMiu +(1-ResultYSchemeII(i))*POperation;
                else
                    synthesisSchemeII(i) = PMiu;
                end

            end

            B = ResultXSchemeII .* (synthesisSchemeII' * ones(1, c));

            B_scheme2Centralized = B;
    
            
            [resultedInterference, exceedInterferenceBound] = checkResultedInference(B_scheme2Centralized, n, m, GtildeAll, infBound);            
        else
                        RetGUROBI_FCC(run) = 0;

        end
        


    
    
        
    
    SchemeIISolutionFlag = 0; % Todo: disable execution of scheme2distributed
        %% Scheme2Distributed
        %-------------------------------|
        %         scheme2distributed    |
        %-------------------------------|
    if(SchemeIISolutionFlag) % only run scheme2distributed when scheme2centralized has solution.

        B = zeros(n, c);
        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
        recordPerf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
        sumUtilityScheme2Distributed = [sumUtility];

        stop = 0;
        Bbackup = B;
        updateCount = 0;
        
        % record the initial SINR for all SUs
        SINRvarianceWhitecat = 0; % sum of variance for all WBS in the whole process of convergence
        delta1step = 0;
        SINR_ETs_previous = [];
        SINRofETs = [];
        [SINRofETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_previous = SINRofETs;
        
        sumUtilityScheme2distributed = [];
        while (stop == 0)
            for i = 1: n
                [B, updateFlag] = update_Scheme2Distributed(seq(i), B, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta, infBound, POperation);
                updateCount = updateCount + updateFlag;
                
                if(updateFlag)  % there is a update
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
                    sumUtilityScheme2distributed(end+1) = sumUtility; % record the trace of sum utility
                    
%                     % calculate the percentage of variance and sum up!
%                     [SINRofETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
%                     if(size(SINRofETs)~= size(SINR_ETs_previous))
%                         dimensionsNoagree = 1;
%                     end
%                     delta1step = sum(abs((SINRofETs - SINR_ETs_previous)./SINR_ETs_previous));
                    
                else
                    sumUtilityScheme2distributed(end+1) = sumUtilityScheme2distributed(end);
                end
                
%                 SINRvarianceScheme2distributed = SINRvarianceScheme2distributed + delta1step;

                if i == n   % calculate utility after dealing with su n's channel, record the utility after one round optimization
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
                    recordPerf = [recordPerf; sumUtility, averageI, averageP, averageSINR, stdSINR];
                end
                SINR_ETs_previous = SINRofETs;
            end
            
            if (updateCount>10)
               needchecking = 1;
            end

            if (isequal(B, Bbackup))	           % B and B_backup(the previous B) are the same!
                stop = 1;
            end
            Bbackup = B; % Bbackup records the current B
        end

%         convergenceStepWhitecat(run) = size(sumUtilityScheme2distributed, 2);
%         SINRvarianceWhitecat_container(run) = SINRvarianceWhitecat;
        
%         dica_perf = recordPerf(end, :);
%         
%         % check sinr on end users.
%         SINR_ETs_Scheme2distributed = []; % there should be n*nET values
%         %SINR_ETs_whitecat = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
%         [SINR_ETs_Scheme2distributed, worstSINR_Scheme2distributed, fair_Scheme2distributed] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);        
%         SINR_ETs_Scheme2distributed_container = [SINR_ETs_Scheme2distributed_container, SINR_ETs_Scheme2distributed];
%         fair_Scheme2distributed_container = [fair_Scheme2distributed_container, fair_Scheme2distributed];
%         worstSINR_Scheme2distributed_container = [worstSINR_Scheme2distributed_container, worstSINR_Scheme2distributed];
%         
% 
%         disp('snrRatio_Scheme2Distributed:');
%         snrRatio_dica = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);    % output quai SINR of all users
        

        B_Scheme2Distributed = B;

    end
else  %% the above codes execute scheme II, both centralized and distributed schemes.


%% start here the channel allocation schemes for ECC (scheme I)


        %%% Initialize channels asignment randomly
        initialB = generateInitialChannelAllocationB_for_schemeI(n,c, P);
                    
        % check sinr on end users.
        SINR_ETs_random = []; % there should be n*nET values
        %SINR_ETs_random = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
        [SINR_ETs_random, worstSINR_random, fair_random] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_random_container = [SINR_ETs_random_container, SINR_ETs_random];
        fair_random_container = [fair_random_container, fair_random];
        worstSINR_random_container = [worstSINR_random_container, worstSINR_random];
        
        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);

        random_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
        disp(':');        
        snrRatio_random = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);    % output quai SINR of all users
        B_random = B;       
        
%         %-----------------------------------------------------------%
%         %         optimation: channel allocation                   %
%         %         - objective function is quadratic                 %
%         %         - constraints are linear                          %
%         %-----------------------------------------------------------%
%         % optimization:
% %                 channelUsage = x(1: n*c);
% %         powerLevelFinal = x(n*c+1: 2*n*c);
%         minPower = 4;
%         [channelUsage, powerLevelFinal] = jointPowerChannelAlloc_GUROBI(n, c, P, Gtilde, delta, minPower);


        %%
        %---------------------------------------------
        %         optimization by GUROBI
        %---------------------------------------------
        % store Gtilde into file
        GtildeInOneRow=[];
        for i=1:n
            GtildeInOneRow = [GtildeInOneRow,Gtilde(i, :)];
        end
        dlmwrite('/Users/max/Documents/git_li/channel-power-allocation-802.22/generated_data/Gtilde.txt', GtildeInOneRow, ' ');

        %---- store powerRatio into file
        a=[];
        condenseP=zeros(n, c);  % P-->condenseP which is n x c: p_i^k, i\in N, k\in C
        for i= 1:n
            condenseP(i,:) = sum(P( (i-1)*c+1 : i*c, :), 2);
            a = [a, condenseP(i,:)]; % finaly get (nc x 1) matrix
        end
        % put two additional constant in the end of the powers
        a=[a, delta]; % add noise 'delta'
        %a=[a, 1/(SUcellRadius^pathlossfactor)];    % signalfade 
        dlmwrite('/Users/max/Documents/git_li/channel-power-allocation-802.22/generated_data/maxPermittedPower.txt', a, ' ');
        %save('/home/li/work/tools/lindo/lindoapi/samples/c/dica/power.txt', a, '-double');
        pause(1);
        
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

        dlmwrite('/Users/max/Documents/git_li/channel-power-allocation-802.22/generated_data/arrayH.txt', arrayH, ' ');
        
        % model.obj
        % N_0/p_{i,k}
        %NoisePowerRatio = delta./(condenseP*SUcellRadius^(-pathlossfactor));
        NoisePowerRatio = delta./(condenseP);
        NoisePowerRatioInOneRow =[];
        for i=1:c
            NoisePowerRatioInOneRow = [NoisePowerRatioInOneRow, NoisePowerRatio(:, i)'];
        end
        dlmwrite('/Users/max/Documents/git_li/channel-power-allocation-802.22/generated_data/NoisePowerRatio.txt', NoisePowerRatioInOneRow, ' ');
        
        % model.A, parameters in constraints
        % n x (n*c)
        % assume n = 4, c = 2, then model.A is,
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
        optimizaionModel1.rhs = ones(1, n);
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

        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
        GUROBI_Perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];

        SINR_ETs_optimization = []; % there should be n*nET values
%         SINR_ETs_optimization = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
        [SINR_ETs_optimization, worstSINR_optimization, fair_optimization] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_optimization_container = [SINR_ETs_optimization_container, SINR_ETs_optimization];
        fair_optimization_container = [fair_optimization_container, fair_optimization];
        worstSINR_optimization_container = [worstSINR_optimization_container, worstSINR_optimization]; 
        disp('snrRatio_optimization');
        snrRatio_optimization = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);     % output quai SINR of all users
        
        
        B_optimization= B;
        
        %%
        %---------------------------|
        %         WhiteCat          |
        %---------------------------|

        B = initialB;
        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
        recordPerf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
        sumUtilityWhitecat = [sumUtility];

        stop = 0;
        Bbackup = B;
        updateCount = 0;
        
        % record the initial SINR for all SUs
        SINRvarianceWhitecat = 0; % sum of variance for all WBS in the whole process of convergence
        delta1step = 0;
        SINR_ETs_previous = [];
        SINRofETs = [];
        [SINRofETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_previous = SINRofETs;
        
        while (stop == 0)
            for i = 1: n
                [B, updateFlag] = update(seq(i), B, P, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta);
                updateCount = updateCount + updateFlag;
                
                if(updateFlag)  % there is a update
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
                    sumUtilityWhitecat(end+1) = sumUtility; % record the trace of sum utility
                    
                    % calculate the percentage of variance and sum up!
                    [SINRofETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
                    if(size(SINRofETs)~= size(SINR_ETs_previous))
                        dimensionsNoagree = 1;
                    end
                    delta1step = sum(abs((SINRofETs - SINR_ETs_previous)./SINR_ETs_previous));
                    
                else
                    sumUtilityWhitecat(end+1) = sumUtilityWhitecat(end);
                end
                
                SINRvarianceWhitecat = SINRvarianceWhitecat + delta1step;

                if i == n   % calculate utility after dealing with su n's channel, record the utility after one round optimization
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
                    recordPerf = [recordPerf; sumUtility, averageI, averageP, averageSINR, stdSINR];
                end
                SINR_ETs_previous = SINRofETs;
            end
            
%             if (updateCount>100)
%                needchecking = 1;
%             end

            if (isequal(B, Bbackup))	           % B and B_backup(the previous B) are the same!
                stop = 1;
            end
            Bbackup = B; % Bbackup records the current B
        end

        convergenceStepWhitecat(run) = size(sumUtilityWhitecat, 2);
        SINRvarianceWhitecat_container(run) = SINRvarianceWhitecat;
        
        dica_perf = recordPerf(end, :);
        
        % check sinr on end users.
        SINR_ETs_whitecat = []; % there should be n*nET values
        %SINR_ETs_whitecat = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
        [SINR_ETs_whitecat, worstSINR_cat, fair_cat] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);        
        SINR_ETs_whitecat_container = [SINR_ETs_whitecat_container, SINR_ETs_whitecat];
        fair_cat_container = [fair_cat_container, fair_cat];
        worstSINR_cat_container = [worstSINR_cat_container, worstSINR_cat];
        

        disp('snrRatio_dica:');
        snrRatio_dica = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);    % output quai SINR of all users
        
        if (max(snrRatio_dica)>50)
           sss=1; 
        end
        B_cat = B;
        

        %%
        %-----------------------------|
        %         WhiteCase           |
        %-----------------------------|
        % sumUtilityWhitecase contains at most 1000*16 records
        B = initialB;

        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
        recordPerf_self = [sumUtility, averageI, averageP, averageSINR, stdSINR];
        sumUtilityWhitecase = [sumUtility];
        
        stop = 0;
        Bbackup = B;
        updateCount = 0;
        deadlockcount = 0;

        % record the initial SINR for all SUs
        SINRvarianceWhitecase = 0; % sum of variance for all WBS in the whole process of convergence
        SINR_ETs_previous = [];
        SINRofETs = [];
        
        [SINRofETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);        
        SINR_ETs_previous = SINRofETs;        
        
        while (stop == 0)
            for i = 1: n
                [B, updateFlag] = selfishUpdate(seq(i), B, P, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta);
                updateCount = updateCount + updateFlag;
                if(updateFlag)
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
                    sumUtilityWhitecase(end+1) = sumUtility; % record the trace of sum utility
                else
                    sumUtilityWhitecase(end+1) = sumUtilityWhitecase(end);
                end

                % calculate the percentage of variance and sum up!
                [SINRofETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
                delta1step = sum(abs((SINRofETs - SINR_ETs_previous)./SINR_ETs_previous));                                    
                
                SINRvarianceWhitecase = SINRvarianceWhitecase + delta1step;
                
                if i == n   % calculate utility after dealing with su n's channel, record the utility after one round optimization
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
                    recordPerf_self = [recordPerf_self;  sumUtility, averageI, averageP, averageSINR, stdSINR];
                end
                SINR_ETs_previous = SINRofETs;                
            end

            if (isequal(B, Bbackup))	           % B and B_backup(the previous B) are the same!
                stop = 1;
            end
            Bbackup = B; % Bbackup records the current B
            if(deadlockcount==50)
                stop =1;
            end 
            deadlockcount = deadlockcount + 1;
        end
                
        convergenceStepWhitecase(run) = size(sumUtilityWhitecase, 2);
        SINRvarianceWhitecase_container(run) = SINRvarianceWhitecase;

        SINR_ETs_whitecase = []; % there should be n*nET values
        % SINR_ETs_whitecase = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s); 
        [SINR_ETs_whitecase, worstSINR_case, fair_case] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_whitecase_container = [SINR_ETs_whitecase_container, SINR_ETs_whitecase];
        fair_case_container = [fair_case_container, fair_case];
        worstSINR_case_container = [worstSINR_case_container, worstSINR_case];        
        
        
        selfishUpdate_perf = recordPerf_self(end, :);
        disp('snrRatio_self');
        snrRatio_self = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);     % output quai SINR of all users
        
        B_case = B;
        
        %------- End of selfish -----------------|

        


        %%        
        %-------------------------------------------------------------------|
        %         noregret based learning algorithm                         |
        %-------------------------------------------------------------------|
        B = initialB;
        sumUtilityNoregret = [];

       
        [B, sumUtilityNoregret, SINRvarianceNoregret] = noregretlearning(seq, B, P, Gtilde, GtildeETsSUs, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta);
        
        convergenceStepNoregret(run) = size(sumUtilityNoregret, 2);
        SINRvarianceNoregret_container(run) = SINRvarianceNoregret;      
        
        SINR_ETs_noregret = []; % there should be n*nET values
        % SINR_ETs_noregret = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
        [SINR_ETs_noregret, worstSINR_noregret, fair_noregret] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);        
        SINR_ETs_noregret_container = [SINR_ETs_noregret_container, SINR_ETs_noregret];
        fair_noregret_container = [fair_noregret_container, fair_noregret];
        worstSINR_noregret_container = [worstSINR_noregret_container, worstSINR_noregret];          
        
        B_noregret=B;
        disp('snrRatio_noregret:');
        snrRatio_noregret = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);     % output quai SINR of all users
        
        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
        noregret_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
        
        B_noregret =B;
        %------- End of noregret learning -----------------|

        %%        
        
        %-------------------------------------------------------------------|
        %         potential game:
        %         minimize the difference between received sigle power and
        %         the produced and received interference
        %         ref: pimrc_2011
        % powerLevels: is the number of power levels.
        %-------------------------------------------------------------------|
        B = initialB;
        powerLevels = 10;
        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
        recordPerf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
        sumUtilityPotentialGame = [sumUtility];

        stop = 0;
        Bbackup = B;
        updateCount = 0;
        
        needchecking=0;
        % record the initial SINR for all SUs
        SINRvariancePotentialGame = 0; % sum of variance for all WBS in the whole process of convergence
        delta1step = 0;
        SINR_ETs_previous = [];
        SINRofETs = [];
        [SINRofETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_previous = SINRofETs;
        
        while (stop == 0)
            for i = 1: n
                [B, updateFlag] = update_distributed_ref(seq(i), B, P, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta, powerLevels);
                updateCount = updateCount + updateFlag;
                
                if(updateFlag)  % there is a update
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
                    sumUtilityPotentialGame(end+1) = sumUtility; % record the trace of sum utility
                    
                    % calculate the percentage of variance and sum up!
                    [SINRofETs, worstSINR, fairness] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
                    if(size(SINRofETs)~= size(SINR_ETs_previous))
                        dimensionsNoagree = 1;
                    end
                    delta1step = sum(abs((SINRofETs - SINR_ETs_previous)./SINR_ETs_previous));
                    
                else
                    sumUtilityPotentialGame(end+1) = sumUtilityPotentialGame(end);
                end
                
                SINRvariancePotentialGame = SINRvariancePotentialGame + delta1step;

                if i == n   % calculate utility after dealing with su n's channel, record the utility after one round optimization
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
                    recordPerf = [recordPerf; sumUtility, averageI, averageP, averageSINR, stdSINR];
                end
                SINR_ETs_previous = SINRofETs;
            end
            
            if (updateCount>1000)
               needchecking = 1;
            end

            if (isequal(B, Bbackup) || needchecking)	           % B and B_backup(the previous B) are the same!
                stop = 1;
            end
            Bbackup = B; % Bbackup records the current B
        end

        convergenceStepPotentialGame(run) = size(sumUtilityPotentialGame, 2);
        SINRvariancePotentialGame_container(run) = SINRvariancePotentialGame;
        
        PotentialGame_perf = recordPerf(end, :);
        B_PotentialGame = B; 
        % check sinr on end users.
        SINR_ETs_PotentialGame = []; % there should be n*nET values
        %SINR_ETs_whitecat = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
        [SINR_ETs_PotentialGame, worstSINR_PotentialGame, fair_PotentialGame] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);       
        SINR_ETs_PotentialGame_container = [SINR_ETs_PotentialGame_container, SINR_ETs_PotentialGame];
        fair_PotentialGame_container = [fair_PotentialGame_container, fair_PotentialGame];
        worstSINR_PotentialGame_container = [worstSINR_PotentialGame_container, worstSINR_PotentialGame];
        

        disp('snrRatio_PotentialGame:');
        snrRatio_PotentialGame = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);    % output quai SINR of all users
        
        if (max(snrRatio_PotentialGame)>50)
           sss=1; 
        end
        B_PotentialGame = B;        
        %---------potential game ends!


        %%
        
        utilityHistory(:, run) = [random_perf(1); dica_perf(1); selfishUpdate_perf(1); noregret_perf(1); GUROBI_Perf(1); PotentialGame_perf(1)];
        powerHistory(:, run) = [random_perf(3); dica_perf(3); selfishUpdate_perf(3); noregret_perf(3); GUROBI_Perf(3); PotentialGame_perf(3)];
        averageSinrHistory(:, run) = [random_perf(4); dica_perf(4); selfishUpdate_perf(4); noregret_perf(4); GUROBI_Perf(4); PotentialGame_perf(4)];
        averageStdHistory(:, run) = [random_perf(5); dica_perf(5); selfishUpdate_perf(5); noregret_perf(5); GUROBI_Perf(5); PotentialGame_perf(5)];
    
        
        
%         %  with lindo
%         utilityHistory(:, run) = [random_perf(1); GUROBI_Perf(1); dica_perf(1); selfishUpdate_perf(1); noregret_perf(1)];
%         powerHistory(:, run) = [random_perf(3); GUROBI_Perf(3); dica_perf(3); selfishUpdate_perf(3); noregret_perf(3)];
%         averageSinrHistory(:, run) = [random_perf(4); GUROBI_Perf(4); dica_perf(4); selfishUpdate_perf(4); noregret_perf(4)];
%         averageStdHistory(:, run) = [random_perf(5); GUROBI_Perf(5); dica_perf(5); selfishUpdate_perf(5); noregret_perf(5)];
        
%         % without Lindo
%         % dica_perf contains: sumUtility, averageI, averageP, averageSINR, stdSINR
%         powerHistory(:, run) = [random_perf(3); dica_perf(3); selfishUpdate_perf(3); noregret_perf(3)];
%         averageSinrHistory(:, run) = [random_perf(4); dica_perf(4); selfishUpdate_perf(4); noregret_perf(4)];
%         averageStdHistory(:, run) = [random_perf(5); dica_perf(5); selfishUpdate_perf(5); noregret_perf(5)];
        
%        convergenplot(sumUtilityWhitecat, sumUtilityWhitecase, sumUtilityNoregret);
end

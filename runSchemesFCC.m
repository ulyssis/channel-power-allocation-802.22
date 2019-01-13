

function [RetGUROBI_FCC, utilityHistoryFCC, powerHistoryFCC, averageSinrHistoryFCC, averageStdHistoryFCC, ...
    SINR_ETs_centralized_FCC_container, SINR_ETs_distributed_FCC_container, fair_centralized_FCC_container, fair_distributed_FCC_container, ...
    worstSINR_centralized_FCC_container, worstSINR_distributed_FCC_container] = ...
    runSchemesFCC(run, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, ...
    utilityHistoryFCC, powerHistoryFCC, averageSinrHistoryFCC, averageStdHistoryFCC, ...
    SINR_ETs_centralized_FCC_container, SINR_ETs_distributed_FCC_container, fair_centralized_FCC_container, fair_distributed_FCC_container, ...
    worstSINR_centralized_FCC_container, worstSINR_distributed_FCC_container,...
    PMiu, POperation, infBound, RetGUROBI_FCC)

seq = randperm(n);
SchemeIISolutionFlag = 0;
B_scheme2Centralized = [];
B_Scheme2Distributed = [];

%% scheme II

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
            filename_B_scheme2Centralized = "/Users/max/Documents/git_li/channel-power-allocation-802.22/B_scheme2Centralized_POperation_" + POperation + '.csv';
            if isfile(filename_B_scheme2Centralized)
                dlmwrite("/Users/max/Documents/git_li/channel-power-allocation-802.22/B_scheme2Centralized_POperation_" + POperation + ".csv", B_scheme2Centralized, '-append');
            else
                dlmwrite("/Users/max/Documents/git_li/channel-power-allocation-802.22/B_scheme2Centralized_POperation_" + POperation + ".csv", B_scheme2Centralized);
            end
            
            [resultedInterference, exceedInterferenceBound] = checkResultedInference(B_scheme2Centralized, n, m, GtildeAll, infBound);            
        else
                        RetGUROBI_FCC(run) = 0;

        end
        

        %% obtain performance
        if (resultStatus)
            % check sinr on end users.
            SINR_ETs_centralized_FCC = []; % there should be n*nET values
            %SINR_ETs_centralized_FCC = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
            [SINR_ETs_centralized_FCC, worstSINR_centralized_FCC, fair_centralized_FCC] = SINR_ETs_cellReSelection(B_scheme2Centralized, n, GtildeETsSUs, nET, TVpower, delta);
            SINR_ETs_centralized_FCC_container = [SINR_ETs_centralized_FCC_container, SINR_ETs_centralized_FCC];
            fair_centralized_FCC_container = [fair_centralized_FCC_container, fair_centralized_FCC];
            worstSINR_centralized_FCC_container = [worstSINR_centralized_FCC_container, worstSINR_centralized_FCC];

            [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B_scheme2Centralized, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);

            centralized_FCC_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
            disp(':');        
            snrRatio_centralized_FCC = output(B_scheme2Centralized, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);    % output quai SINR of all users

        end
    
    
        
    
        %% Scheme2Distributed
        %-------------------------------|
        %         scheme2distributed    |
        %-------------------------------|
    if(SchemeIISolutionFlag) % only run scheme2distributed when scheme2centralized has solution.

        B = zeros(n, c);
        recordPerf = [];

        stop = 0;
        Bbackup = B;
        updateCount = 0;
        deadlockcount = 0;
        
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
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);
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
                    [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);
                    recordPerf = [recordPerf; sumUtility, averageI, averageP, averageSINR, stdSINR];
                end
                SINR_ETs_previous = SINRofETs;
            end
            
            deadlockcount = deadlockcount + 1;
            if (updateCount>10 )
               needchecking = 1;
            end

            if (isequal(B, Bbackup) || deadlockcount > 50 )	           % B and B_backup(the previous B) are the same!
                stop = 1;
            end
            Bbackup = B; % Bbackup records the current B
        end
        
        dlmwrite("/Users/max/Documents/git_li/channel-power-allocation-802.22/sumUtilityScheme2distributed.csv", sumUtilityScheme2distributed, '-append');


        B_Scheme2Distributed = B;
        
        filename_B_scheme2Distributed = "/Users/max/Documents/git_li/channel-power-allocation-802.22/B_scheme2Distributed_POperation_" + POperation + '.csv';
        if isfile(filename_B_scheme2Distributed)
            dlmwrite("/Users/max/Documents/git_li/channel-power-allocation-802.22/B_scheme2Distributed_POperation_" + POperation + ".csv", B_Scheme2Distributed, '-append');
        else
            dlmwrite("/Users/max/Documents/git_li/channel-power-allocation-802.22/B_scheme2Distributed_POperation_" + POperation + ".csv", B_Scheme2Distributed);
        end
        
        %% obtain performance of distributed scheme
        % check sinr on end users.
        SINR_ETs_distributed_FCC = []; % there should be n*nET values
        %SINR_ETs_distributed_FCC = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
        [SINR_ETs_distributed_FCC, worstSINR_distributed_FCC, fair_distributed_FCC] = SINR_ETs_cellReSelection(B_Scheme2Distributed, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_distributed_FCC_container = [SINR_ETs_distributed_FCC_container, SINR_ETs_distributed_FCC];
        fair_distributed_FCC_container = [fair_distributed_FCC_container, fair_distributed_FCC];
        worstSINR_distributed_FCC_container = [worstSINR_distributed_FCC_container, worstSINR_distributed_FCC];
        
        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B_Scheme2Distributed, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);

        distributed_FCC_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
        disp(':');        
        snrRatio_distributed_FCC = output(B_Scheme2Distributed, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);    % output quai SINR of all users

        %% record performace of the two schemes of scheme II
        utilityHistoryFCC(:, run) = [centralized_FCC_perf(1); distributed_FCC_perf(1)];
        powerHistoryFCC(:, run) = [centralized_FCC_perf(3); distributed_FCC_perf(3)];
        averageSinrHistoryFCC(:, run) = [centralized_FCC_perf(4); distributed_FCC_perf(4)];
        averageStdHistoryFCC(:, run) = [centralized_FCC_perf(5); distributed_FCC_perf(5)];

    end



        




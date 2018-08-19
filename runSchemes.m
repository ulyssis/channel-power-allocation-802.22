       
function [utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_lindo_container, SINR_ETs_lindo2_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_lindo_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_lindo_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, B_random, B_cat, B_case, B_lindo, B_noregret, B_PotentialGame, B_lindoCAPA, flag_resolve] = runSchemes(run, P, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_lindo_container, SINR_ETs_lindo2_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_lindo_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_lindo_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container)
%function [powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_lindo_container, SINR_ETs_lindo2_container, SINR_ETs_noregret_container, fair_random_container, fair_cat_container, fair_case_container, fair_noregret_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, snrRatio_random, snrRatio_dica, snrRatio_self, snrRatio_noregret] = runSchemes(run, P, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_lindo_container, SINR_ETs_lindo2_container, SINR_ETs_noregret_container, fair_random_container, fair_cat_container, fair_case_container, fair_noregret_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container)


seq = randperm(n);
   
        % Initialize channels asignment randomly
        B = zeros(n, c);
        doagain=1;
            while (doagain)
                for i = 1 : n
                   B(i, :) = P((i-1)*c + floor((1+c * rand)), :);        
                end
                doagain=0;
                for i=1:c
                   if (nnz(B(:,i))==1)
                       doagain=1;
                       break;
                   end
                end
            end

        % check sinr on end users.
        SINR_ETs_random = []; % there should be n*nET values
        %SINR_ETs_random = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
        [SINR_ETs_random, worstSINR_random, fair_random] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_random_container = [SINR_ETs_random_container, SINR_ETs_random];
        fair_random_container = [fair_random_container, fair_random];
        worstSINR_random_container = [worstSINR_random_container, worstSINR_random];
        
        
        initialB = B;       % record the initial B for selfishUpdate
        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);

        random_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
        disp(':');        
        snrRatio_random = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);    % output quai SINR of all users
        B_random = B;

        
        
%         %---------------------------|
%         %         optimation: joint power and channel allocation
%         %         objective function is quadratic
%         %         constraints are linear
%         %---------------------------|
%         % optimization:
% %                 channelUsage = x(1: n*c);
% %         powerLevelFinal = x(n*c+1: 2*n*c);
%         minPower = 4;
%         [channelUsage, powerLevelFinal] = jointPowerChannelAlloc_gurobi(n, c, P, Gtilde, delta, minPower);


        %%        
        %---------------------------------------------
        %         parameters input for OPT solver
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
        
       
        
        % generate compounded num_constants for lindoapi
        % h_{i,j}*p_{j,k}/p_{i,k}
        compoundMatrix = [];
        compoundMatrixInOneRow = [];
        % the compoundMatrix is a nn x c matrix.
        % 
        % p_i^k, i\in N, k\in C
        % ./
        % p_1^k, k\in C 
        % .* 
        % h_{1, i}, i\in N                 nxc
        % ----------------------
        % p_i^k, i\in N, k\in C
        % ./
        % p_2^k, k\in C
        % .* 
        % h_{2, i}, i\in N                 nxc
        % ----------------------
        % .
        % .
        % .
        % ----------------------
        % p_i^k, i\in N, k\in C
        % ./
        % p_n^k, k\in C
        % .* 
        % h_{n, i}, i\in N                 nxc
  
        for i=1:n
            newpart = condenseP./(ones(n, 1)*condenseP(i, :)); % a complete power matrix for all nodes ./ one powermatrix of one node
            newpart = newpart.*(Gtilde(i, :)'*ones(1, c)); % multiply the passloss between any pair.            
            compoundMatrix = [compoundMatrix; newpart];
        end
        for i=1:n*n
            compoundMatrixInOneRow = [compoundMatrixInOneRow, compoundMatrix(i, :)]; 
        end
        dlmwrite('/Users/max/Documents/git_li/channel-power-allocation-802.22/generated_data/compoundCoefficient.txt', compoundMatrixInOneRow, ' ');
        
        % generate compounded num_constants for lindoapi
        % N_0/p_{i,k}
        %NoisePowerRatio = delta./(condenseP*SUcellRadius^(-pathlossfactor));
        NoisePowerRatio = delta./(condenseP);
        NoisePowerRatioInOneRow =[];
        for i=1:n
            NoisePowerRatioInOneRow = [NoisePowerRatioInOneRow, NoisePowerRatio(i, :)];
        end
        dlmwrite('/Users/max/Documents/git_li/channel-power-allocation-802.22/generated_data/NoisePowerRatio.txt', NoisePowerRatioInOneRow, ' ');
        
% % %         %--------The input parameters for Lindo is ready --------%
% % %         % joint power and channel allocation, with lindo
% % %          readMatrixB_jointPowerChannelAllocation(); % call lindo in .c file, to run a script to intriger lindo
% % %          pause(3); % wait for 1s for the result to be written.
% % % 
% % % % if load('resolve') ==1, which indicates lindo returns 'good' results for joint channel-power allocation          
% % % flag_resolve = load('flag_resolve');
% % % 
% % % 
% % % matrixB_vriablePower = load('matrixB_vriablePower');
% % % matrixP_vriablePower = load('matrixP_vriablePower');
% % %          B_lindoCAPA= matrixB_vriablePower.*matrixP_vriablePower;
B_lindoCAPA = [];
        
        
%         % % % % % % % % % % % % % % % % % % % % % % %
%         %   Partial Optimality
%         %   stop here, input Q from Lindo !
%         % % % % % % % % % % % % % % % % % % % % % % %
%         B = condenseP.*Q;
%         Blindo = B;
%         [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
%         lindo_Perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
% 
%         SINR_ETs_lindo = []; % there should be n*nET values
% %         SINR_ETs_lindo = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
%         [SINR_ETs_lindo, worstSINR_lindo, fair_lindo] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
%         SINR_ETs_lindo_container = [SINR_ETs_lindo_container, SINR_ETs_lindo];
%         disp('snrRatio_lindo');
%         snrRatio_lindo = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);     % output quai SINR of all users

       
        
        % % % % % % % % % % % % % % % % % % % % % % %
        %   Global Optimality
        %   stop here, input B from Lindo !
        % % % % % % % % % % % % % % % % % % % % % % %
        readMatrixB(); % call lindo in .c file, to run a script to intriger lindo
        B= condenseP.*load('matrixB');
        [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
        lindo_Perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];

        SINR_ETs_lindo = []; % there should be n*nET values
%         SINR_ETs_lindo = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s);
        [SINR_ETs_lindo, worstSINR_lindo, fair_lindo] = SINR_ETs_cellReSelection(B, n, GtildeETsSUs, nET, TVpower, delta);
        SINR_ETs_lindo_container = [SINR_ETs_lindo_container, SINR_ETs_lindo];
        fair_lindo_container = [fair_lindo_container, fair_lindo];
        worstSINR_lindo_container = [worstSINR_lindo_container, worstSINR_lindo]; 
        disp('snrRatio_lindo');
        snrRatio_lindo = output(B, Gtilde, GtildeAll, n, m, TVpower, SUcellRadius, delta, pathlossfactor);     % output quai SINR of all users
        
        
        B_lindo= B;
        



    
% %        call lindo with c language, get B!
% %           mex autoGenerated.c -o autoGenerated -L/home/li/work/tools/lindo/lindoapi/bin/linux64 -llindo64  -lmosek64 -llindojni -lconsub3 -lc -ldl -lm -lguide -lpthread -lsvml -limf -lirc
% %         ------- End of lindo for joint power/channel allocation -----------------|
        
        
        
        
%%
        %---------------------------|
        %         WhiteCat          |
        %---------------------------|

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
        B_dica = B; 
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
        B_selfish = B; 
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
        
        
        % with lindo(fixed power), lindo(variable power), and distributed
        % variable power
        utilityHistory(:, run) = [random_perf(1); dica_perf(1); selfishUpdate_perf(1); noregret_perf(1); lindo_Perf(1)];
        powerHistory(:, run) = [random_perf(3); dica_perf(3); selfishUpdate_perf(3); noregret_perf(3); lindo_Perf(3)];
        averageSinrHistory(:, run) = [random_perf(4); dica_perf(4); selfishUpdate_perf(4); noregret_perf(4); lindo_Perf(4)];
        averageStdHistory(:, run) = [random_perf(5); dica_perf(5); selfishUpdate_perf(5); noregret_perf(5); lindo_Perf(5)];
        
        
        
%         %  with lindo
%         utilityHistory(:, run) = [random_perf(1); lindo_Perf(1); dica_perf(1); selfishUpdate_perf(1); noregret_perf(1)];
%         powerHistory(:, run) = [random_perf(3); lindo_Perf(3); dica_perf(3); selfishUpdate_perf(3); noregret_perf(3)];
%         averageSinrHistory(:, run) = [random_perf(4); lindo_Perf(4); dica_perf(4); selfishUpdate_perf(4); noregret_perf(4)];
%         averageStdHistory(:, run) = [random_perf(5); lindo_Perf(5); dica_perf(5); selfishUpdate_perf(5); noregret_perf(5)];
        
%         % without Lindo
%         % dica_perf contains: sumUtility, averageI, averageP, averageSINR, stdSINR
%         powerHistory(:, run) = [random_perf(3); dica_perf(3); selfishUpdate_perf(3); noregret_perf(3)];
%         averageSinrHistory(:, run) = [random_perf(4); dica_perf(4); selfishUpdate_perf(4); noregret_perf(4)];
%         averageStdHistory(:, run) = [random_perf(5); dica_perf(5); selfishUpdate_perf(5); noregret_perf(5)];
        
%        convergenplot(sumUtilityWhitecat, sumUtilityWhitecase, sumUtilityNoregret);


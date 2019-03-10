function [utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
    SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, ...
    B_random, B_cat, B_case, B_optimization, B_noregret, B_PotentialGame] ...
    = runSchemes(run, w, P, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory, powerHistory, ...
    averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
    SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, PMiu)


% matrix B is expanded for w times
% now B's size is n*w X c

seq = randperm(n);
   %% random channel allocation
   
        % Initialize channels asignment randomly
        B = zeros(n*w, c);
        doagain=1;
            while (doagain)
                for i = 1 : n*w
                   B(i, :) = P(floor((1+c * rand)), :);        
                end
                doagain=0;
                for i=1:c
                   if (nnz(B(:,i))==1)
                       doagain=1;
                       break;
                   end
                end
            end

       initialB = B;       % record the initial B for selfishUpdate



        
        %---------------------------------------------
        %         optimation: channel allocation                   
        %         - objective function is quadratic                 
        %         - constraints are linear                         
        %         - optimization by GUROBI
        %---------------------------------------------

        B = GUROBI_ECC(n, c, w, P, Gtilde, delta);

        averageP = mean(sum(B, 2));
        [averageShannonCPerCell] = capacityOnETs(B, n, w, GtildeETsSUs, nET);
        B_optimization = B;
        
        %%
        %---------------------------|
        %         WhiteCat          |
        %---------------------------|

        B = initialB;
%         [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(w, B, n*w, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);
%         recordPerf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
%         sumUtilityWhitecat = [sumUtility];

        stop = 0;
        Bbackup = B;
        updateCount = 0;
        
        while (stop == 0)
            for i = 1: n
                [B, updateFlag] = update(seq(i), w, B, P, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta);
                updateCount = updateCount + updateFlag;
                
%                 if(updateFlag)  % there is a update
%                     [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(w, B, n*w, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);
%                     sumUtilityWhitecat(end+1) = sumUtility; % record the trace of sum utility
%        
%                 else
%                     sumUtilityWhitecat(end+1) = sumUtilityWhitecat(end);
%                 end
                


%                 if i == n   % calculate utility after dealing with su n's channel, record the utility after one round optimization
%                     [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(w, B, n*w, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, PMiu);
%                     recordPerf = [recordPerf; sumUtility, averageI, averageP, averageSINR, stdSINR];
%                 end

            end
            
%             if (updateCount>100)
%                needchecking = 1;
%             end

            if (isequal(B, Bbackup))	           % B and B_backup(the previous B) are the same!
                stop = 1;
            end
            Bbackup = B; % Bbackup records the current B
        end

        %convergenceStepWhitecat(run) = size(sumUtilityWhitecat, 2);
        
        %dica_perf = recordPerf(end, :);

        
%         if (max(snrRatio_dica)>50)
%            sss=1; 
%         end
        condenseB = condense(B, w);
        averageP = mean(sum(condenseB, 2));
        [averageShannonCPerCell] = capacityOnETs(condenseB, n, w, GtildeETsSUs, nET);
        B_cat = condenseB;
        

        %%
        %-----------------------------|
        %         WhiteCase           |
        %-----------------------------|

        
        %------- End of selfish -----------------|

        


        %%        
        %-------------------------------------------------------------------|
        %         noregret based learning algorithm                         |
        %-------------------------------------------------------------------|
   
        %------- End of noregret learning -----------------|

        %%        
        
        %-------------------------------------------------------------------|
        %         potential game:
        %         minimize the difference between received sigle power and
        %         the produced and received interference
        %         ref: pimrc_2011
        % powerLevels: is the number of power levels.
        %-------------------------------------------------------------------|
         
        %---------potential game ends!

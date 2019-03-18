function [centralized_TxPower_allWBSs_allRuns, decentralized_TxPower_allWBSs_allRuns, ...
    centralized_CellThrought_allWBSs_allRuns, decentralized_CellThrought_allWBSs_allRuns, ...
    random_TxPower_allWBSs_allRuns, random_CellThrought_allWBSs_allRuns] ...
    = runSchemes(run, w, P_CVX, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, ...
    SUcellRadius, delta, pathlossfactor, eta, PMiu, ...
    decentralized_TxPower_allWBSs_allRuns, centralized_TxPower_allWBSs_allRuns, ...
    decentralized_CellThrought_allWBSs_allRuns, centralized_CellThrought_allWBSs_allRuns, ...
    random_TxPower_allWBSs_allRuns, random_CellThrought_allWBSs_allRuns)


% matrix B is expanded for w times
% now B's size is n*w X c

seq = randperm(n*w);
   %% random channel allocation
   
   condenseP_CVX = condenseP(P_CVX, n, c);
        % Initialize channels asignment randomly
        B = zeros(n, c);
        doagain=1;

            while (doagain)

                for i = 1 : n
                    %randomly delete c-w colum in Vindex
                    Vindex = 1:c;
                    deletedChannels =  randsample(Vindex, c-w, false);
                    Vindex(deletedChannels) = 0;

                    B(i, :) = condenseP_CVX(i, :).*Vindex;        
                end

                doagain=0;
                
                for i=1:c
                   if (nnz(B(:,i))==0)
                       doagain=1;
                       break;
                   end
                end
            end
            

       initialB = expandB(B, w);       % record the initial B for selfishUpdate


        %---------------------------|
        %         random            |
        %---------------------------|

        randomB = condense(initialB, w);
        [averageShannonCPerCell] = capacityOnETs(randomB, n, w, GtildeETsSUs, nET, delta);
        random_TxPower_allWBSs_allRuns(run, :) = sum(randomB, 2)';
        random_CellThrought_allWBSs_allRuns(run, :) = averageShannonCPerCell;
        
        %---------------------------------------------
        %         optimation: channel allocation                   
        %         - objective function is quadratic                 
        %         - constraints are linear                         
        %         - optimization by GUROBI
        %---------------------------------------------

        B = GUROBI_ECC(n, c, w, P_CVX, Gtilde, delta);

        pause(1); % to avoid failure in capacityOnETs

        disp(B);
        tol = 1.e-6;
        B(B<0 & B> -tol) = 0;
        [averageShannonCPerCell] = capacityOnETs(B, n, w, GtildeETsSUs, nET, delta);
        %B_optimization = B;
        centralized_TxPower_allWBSs_allRuns(run, :) = sum(B, 2)';
        centralized_CellThrought_allWBSs_allRuns(run, :) = averageShannonCPerCell;
        
        B_centralized = B;
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


        numUpdatedWBSs = 0;
                    numLoop = 0;
        while (stop == 0)


            for i = 1: n
                [B, updateFlag] = update(seq(i), w, B, P_CVX, Gtilde, m, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor, eta);
                numUpdatedWBSs = numUpdatedWBSs + updateFlag;
                
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
            numLoop = numLoop + 1;
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
        averageP = sum(condenseB, 2)';
        [averageShannonCPerCell] = capacityOnETs(condenseB, n, w, GtildeETsSUs, nET, delta);
        %B_cat = condenseB;
        decentralized_TxPower_allWBSs_allRuns(run, :) = averageP;
        decentralized_CellThrought_allWBSs_allRuns(run, :) = averageShannonCPerCell;

        B_decentralzied = condenseB;
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

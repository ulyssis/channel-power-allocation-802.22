%---------------------------------------------------------
%-     Channel seclection in cognitive radio networks.
%---------------------------------------------------------
% distributed greedy search, feedback infomation is used in the utility
% each user looks for the channel to minimize its utility.
% 4 schemes: 'random', 'update' and 'selfishUpdate', 'noregret'. 
% Di Li & James Gross, Di Li, 28/10/11
% Di Li,               Di Li, 28/8/18


% utility fucntion: U = f_i/P_i + \sum (f_{ij}/P_j), where f is the
% interferece, P_i is the transmission power, Channel c is ommited here for clearity.

% run distributed schemes and centralized scheme with GUROBI

% to use cvx:
% cd /Users/max/Documents/git_li/channel-power-allocation-802.22/cvx
% cvx_setup

% to use GUROBI
% cd /Library/gurobi801/mac64/matlab
% gurobi_setup

% -- old
% if lindo is called, the following step is needed to automize the running.
% cd /home/li/work/dev/DiCAPS/2011OtcGameOptReform/codes_cameraready
% mex readMatrixB.cpp 
% mex readMatrixB_jointPowerChannelAllocation.cpp
% -- old

% name of save fig:
% ECC_n_c_infBound_SUcellRadius_maxP_miniP_runtimes
% FCC_n_c_infBound_SUcellRadius_POperation_maxP_miniP_runtimes

close all;
echo off;
clc;
addpath("/Users/max/Documents/git_li/channel-power-allocation-802.22/cvx");
addpath("/Library/gurobi801/mac64/matlab");
gurobi_setup;
savepath;

    runtimes =  5;  % number of simulation run
    n = 16;    % number of WBS
    %c = 5;     % number of channels, remeber to modify cvx_statusMsg whose length should be c. 
    %m = c;     % number of primary users, with the same number of channels 
    delta = 1*10.^(-13);   % Noise;untitled.eps
    lengthSide = 60000;
    infBound = 5*10.^(-8);     % The interfernce threshold on PU contour  
    % working parameters:
    % ECC, ie-7
    % FCC, 1e-8
    TVpower = 0;
    pathlossfactor = 2;    
    miniP = 1; % xx dbm, the minmum power for users, 30dbm - 1W
    maxP = 10; % xx dbm
    nET = 10;  % number of endterminals in each WBS 
    s = 8; % set standard deviation
    coverage = lengthSide/4/2 * 0.45; % the maximal distance away from the WBS, whihc a terminal can have 
                                     % This value should consider SUcellRadius.
    eta= 1; % the discount of the sum of interference from different WBSs, to represent the interference on the measurement point 

    runSchemesForECC = 0;
    tic;
    
averagePowerOverNumOfChannels = [];
averagePowerCIOverNumOfChannels = [];
averageETSINROverNumOfChannels = [];
averageETSINRCIOverNumOfChannels = [];

for c = 5:1:5   
    m = c;
    for SUcellRadius = 3000:1000:3000 % 1000:1000:7000


    utilityHistory=[];
    powerHistory=[];
    averageSinrHistory = [];
    averageStdHistory = [];

    SINR_ETs_random_container = [];     
    SINR_ETs_whitecat_container = [];
    SINR_ETs_whitecase_container = [];    
    SINR_ETs_noregret_container = [];  
    SINR_ETs_PotentialGame_container = [];  
    SINR_ETs_optimization_container = [];
    SINR_ETs_optimization2_container = [];

    fair_random_container = [];     
    fair_cat_container = [];     
    fair_case_container = [];     
    fair_optimization_container =[];
    fair_noregret_container = [];  
    fair_PotentialGame_container = [];  

    worstSINR_random_container = [];
    worstSINR_cat_container = [];     
    worstSINR_case_container = [];     
    worstSINR_optimization_container = [];     
    worstSINR_noregret_container = [];   
    worstSINR_PotentialGame_container = []; 

    NOperatingWBSs = [];
    utilityHistoryFCC = [];
    powerHistoryFCC = [];
    averageSinrHistoryFCC = [];
    averageStdHistoryFCC = [];
    SINR_ETs_centralized_FCC_container = []; 
    SINR_ETs_distributed_FCC_container = [];
    fair_centralized_FCC_container = []; 
    fair_distributed_FCC_container = [];
    worstSINR_centralized_FCC_container = [];
    worstSINR_distributed_FCC_container = [];

    convergenceStepWhitecat = zeros(1, runtimes);
    convergenceStepWhitecase = zeros(1, runtimes);
    convergenceStepNoregret = zeros(1, runtimes);
    convergenceStepPotentialGame = zeros(1, runtimes);

    SINRvarianceWhitecat_container = zeros(1, runtimes);
    SINRvarianceWhitecase_container = zeros(1, runtimes);
    SINRvarianceNoregret_container = zeros(1, runtimes);
    SINRvariancePotentialGame_container = zeros(1, runtimes);

    B_random =[];
    B_cat =[];
    B_case=[];
    B_optimization=[];
    B_noregret=[];

    lp_container=[];
    cvx_container=[];

    RetGUROBI_FCC = zeros(1, runtimes); % record whether the execution of GUROBI for scheme II obtains solution.

        baseDir = '/Users/max/Documents/git_li/channel-power-allocation-802.22/';
        FileNameSumUtilityScheme2distributed = fullfile(baseDir, 'sumUtilityScheme2distributed.csv');
        delete(FileNameSumUtilityScheme2distributed);

    for POperation = 4:1:4
        PMiu = 1/POperation;

        contains = dir(baseDir);
        for k = 1:length(contains)
            if (strcmp(contains(k).name, char("B_scheme2Centralized_POperation_" + POperation + ".csv")))
                fullFileName = fullfile(baseDir, contains(k).name);
                delete(fullFileName);
            end
        end



    for run = 1: runtimes % the number of simulations


    %% plot the distibution of stations and nodes
    %    plotlocation(n, m, lengthSide, posSU, posET, posTVContor);

        if(runSchemesForECC)

        %%---------- cvx ---------------
        %    decide the maximal transmission power by solving the convex problem with cvx 
        %    while loop is used to generate solution feasible topologies
            max_cvx_statusMsg =2;
            while(max_cvx_statusMsg > 1) % contains cvx_optval which is the value of the objective function 
                [posSU, posET, posTVContor, Gtilde, GtildeETsSUs, GtildeAll] = geoinfo(n, m, nET, lengthSide, coverage, SUcellRadius, pathlossfactor, s);            
                [P_CVX, max_cvx_statusMsg] = maximalPowerPlanningCVX(n, m, infBound, GtildeAll, miniP, maxP);
            end
                cvx=sum(P_CVX,2);
                cvx_container = [cvx_container, cvx];

                %% run channel assignment scheme I and comparison schemes, which are designed for FCC
                [utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
                    SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, ...
                    B_random, B_cat, B_case, B_optimization, B_noregret, B_PotentialGame] ...
                    = runSchemes(run, P_CVX, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory, powerHistory, ...
                    averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
                    SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, PMiu);
        else
        %% run channel assignment scheme II for FCC
         [RetGUROBI_FCC, utilityHistoryFCC, powerHistoryFCC, averageSinrHistoryFCC, averageStdHistoryFCC, NOperatingWBSs, ...
            SINR_ETs_centralized_FCC_container, SINR_ETs_distributed_FCC_container, fair_centralized_FCC_container, fair_distributed_FCC_container, ...
            worstSINR_centralized_FCC_container, worstSINR_distributed_FCC_container] = ...
            runSchemesFCC(run, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, ...
            utilityHistoryFCC, powerHistoryFCC, averageSinrHistoryFCC, averageStdHistoryFCC, NOperatingWBSs, ...
            SINR_ETs_centralized_FCC_container, SINR_ETs_distributed_FCC_container, fair_centralized_FCC_container, fair_distributed_FCC_container, ...
            worstSINR_centralized_FCC_container, worstSINR_distributed_FCC_container,...
            PMiu, POperation, infBound, RetGUROBI_FCC);
        end

    end

    %plotLog = POperation;
    plotLog = SUcellRadius/100;

        if(runSchemesForECC)
        %      plots of scheme I
                printplotsCAschemes1(plotLog, n, nET, ...
                    powerHistory, averageSinrHistory, ...
                    SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
                    SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container); 

        else
        %      plots of schemII, both centralized and distributed
                printplotsCAschemes2(plotLog, n, nET, ...
                    utilityHistoryFCC, powerHistoryFCC, averageSinrHistoryFCC, NOperatingWBSs, ...
                    SINR_ETs_centralized_FCC_container, SINR_ETs_distributed_FCC_container); 

        %% sumUtility of schemeII distributed.
        %         sumUtilityScheme2distributedAllRuns = load('/Users/max/Documents/git_li/channel-power-allocation-802.22/sumUtilityScheme2distributed.csv');
        %         sumUtilityScheme2distributed = sumUtilityScheme2distributedAllRuns(end, :);
        %         figure(plotLog + 6);
        %         handle1 = plot(sumUtilityScheme2distributed);
        %         set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');

        end
    end      

    %         printplotsWorst20_CAschemes2power(plotLog, n, worstSINR_random_container,worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, worstSINR_optimization_container);

    %         % plot figures for only maximal power decision
    %         worst20(n, worstSINR_random_container,worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, worstSINR_optimization_container);
    %         worst20(n, worstSINR_random_container2,worstSINR_cat_container2, worstSINR_case_container2, worstSINR_noregret_container2, worstSINR_optimization_container2);
    end

    if(runSchemesForECC)
        % averageDataOverNumOfChannels:
        averagePowerOverNumOfChannels = [averagePowerOverNumOfChannels, mean(powerHistory,2)];
        averagePowerCIOverNumOfChannels = [averagePowerCIOverNumOfChannels, 1.96*std(powerHistory,1,2)/sqrt(n)];

        averageETSINROverNumOfChannels = [averageETSINROverNumOfChannels, mean(averageSinrHistory,2)];
        averageETSINRCIOverNumOfChannels = [averageETSINRCIOverNumOfChannels, 1.96*std(averageSinrHistory,1,2)/sqrt(n*nET)];
    end
end


%% plot averageDataOverNumOfChannels and averageDataOverNumOfChannels
if(runSchemesForECC)
    % Average Transmisson Power
    figure(c*10 + 6);
    h = gobjects(size(averagePowerOverNumOfChannels, 1),1);
    for i = 1: size(averagePowerOverNumOfChannels, 1)
        x = [2+0.03*i: 1 :c+0.03*i];
        y = averagePowerOverNumOfChannels(i, :);
        err = averagePowerCIOverNumOfChannels(i, :);
        h(i) = errorbar(x,y,err);
        hold on;
    end
    legend(h, {'Optimization', 'Random Allocation', 'Potential Game', 'No-Regret Learning', 'WhiteCase','whiteCat'}, 'Location','southwest', 'FontSize', 12, 'Color', 'w', 'Box', 'on', 'EdgeColor', 'none');
    xticks(2:1:5);
    xlabel('Number of Available Channels');
    ylabel('Average Tx Power');
    applyhatch(gcf,'|-+.\/');


    % Average SINR
    figure(c*10 + 7);
    h = gobjects(size(averageETSINROverNumOfChannels, 1),1);
    for i = 1: size(averageETSINROverNumOfChannels, 1)
        x = [2+0.03*i: 1 :c+0.03*i];
        y = averageETSINROverNumOfChannels(i, :);
        err = averageETSINRCIOverNumOfChannels(i, :);
        h(i) = errorbar(x,y,err);
        hold on;
    end
    legend(h, {'Optimization', 'Random Allocation', 'Potential Game', 'No-Regret Learning', 'WhiteCase','whiteCat'}, 'Location','northwest', 'FontSize', 12, 'Color', 'w', 'Box', 'on', 'EdgeColor', 'none');
    xticks(2:1:5);
    xlabel('Number of available channels');
    ylabel('Average SINR on End Terminals');
    applyhatch(gcf,'|-+.\/');
end
% %         % Record the sum of utility in the converging
% %         % process in one run
% %         % how to use it: set a breakpoint in the end of function
% %         % "runSchemes", then input the following function in the command
% %         % window
% %         convergenplot(sumUtilityWhitecat, sumUtilityWhitecase, sumUtilityNoregret);
% %%%%%----4 schemes section------------
        
        
% % % % %this is undone
% % % % plotConvergenceSteps(convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret);
% % % % plotSINRVariance(SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container);    
        
% % % % % % 95% confidence interval:
% % % % % % std(convergenceStepNoregret)/sqrt(100)*2*1.962
    toc;

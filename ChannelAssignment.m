%---------------------------------------------------------
%-     Channel seclection in cognitive radio networks.
%---------------------------------------------------------
% distributed greedy search, feedback infomation is used in the utility
% each user looks for the channel to minimize its utility.
% 4 schemes: 'random', 'update' and 'selfishUpdate', 'noregret'. 
% Di Li & James Gross,              28/10/11
% Di Li FCC Schemes                 28/8/18
% Di Li multiple channel            05/03/19

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

    runtimes =  30;  % number of simulation run
    n = 9;    % number of WBS
    c = 4;
    maxNumMultiChannel = 4;     % maximal number of channels allowed to be used 
    minNumMultiChannel = 1;     % minimal number of channels allowed to be used 
    delta = 1*10.^(-13);   % Noise;untitled.eps
    lengthSide = 60000;

    infBound = 1*10.^(-8);         % n=16->  infBound=5*10.^(-8);
                                    % n=9->  infBound=1*10.^(-8);
                                    % The interfernce threshold on PU contour  
    
    TVpower = 0;
    pathlossfactor = 2;    
    miniP = 1; % xx dbm, the minmum power for users, 30dbm - 1W
    maxP = 4; % xx dbm
    nET = 10;  % number of endterminals in each WBS 
    s = 8; % set standard deviation
    coverage = lengthSide/4/2 * 0.5; % the maximal distance away from the WBS, whihc a terminal can have 
                                     % This value should consider SUcellRadius.
    eta= 1; % the discount of the sum of interference from different WBSs, to represent the interference on the measurement point 

    runSchemesForECC = 1;
    
    tic;
    
averagePowerOverNumOfChannels = [];
averagePowerCIOverNumOfChannels = [];
averageETSINROverNumOfChannels = [];
averageETSINRCIOverNumOfChannels = [];

    % Performance in each cell, over all runs.
    decentralized_TxPower_allWBSs_allRuns = zeros(runtimes, n); 
    centralized_TxPower_allWBSs_allRuns = zeros(runtimes, n);
    dyspan14_TxPower_allWBSs_allRuns  = zeros(runtimes, n);
    decentralized_CellThrought_allWBSs_allRuns = zeros(runtimes, n);
    centralized_CellThrought_allWBSs_allRuns = zeros(runtimes, n);
    dyspan14_CellThrought_allWBSs_allRuns  = zeros(runtimes, n);

    m = c;
    for SUcellRadius = 1000:1:1000 % 1000:1000:7000

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
        
    % return a array of empty matrix
    simSesultCell = cell(maxNumMultiChannel, 6); 
    xstick = minNumMultiChannel: 1: maxNumMultiChannel;
for w = xstick
    for run = 1: runtimes % the number of simulations




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
    
                % plot the distibution of stations and nodes
            plotlocation(n, m, lengthSide, posSU, posET, posTVContor);
 
            [centralized_TxPower_allWBSs_allRuns, decentralized_TxPower_allWBSs_allRuns, ...
                dyspan14_TxPower_allWBSs_allRuns, ...
                centralized_CellThrought_allWBSs_allRuns, decentralized_CellThrought_allWBSs_allRuns, ...
                dyspan14_CellThrought_allWBSs_allRuns] ...
                = runSchemes(run, w, maxNumMultiChannel, P_CVX, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, ...
                                SUcellRadius, delta, pathlossfactor, eta, PMiu,...
                                centralized_TxPower_allWBSs_allRuns, decentralized_TxPower_allWBSs_allRuns, ...
                dyspan14_TxPower_allWBSs_allRuns, ...
                centralized_CellThrought_allWBSs_allRuns, decentralized_CellThrought_allWBSs_allRuns,...
                dyspan14_CellThrought_allWBSs_allRuns);      
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
    
            simSesultCell(w, :) = {centralized_TxPower_allWBSs_allRuns,decentralized_TxPower_allWBSs_allRuns, ...
                centralized_CellThrought_allWBSs_allRuns, decentralized_CellThrought_allWBSs_allRuns, ...
                dyspan14_TxPower_allWBSs_allRuns, dyspan14_CellThrought_allWBSs_allRuns};
             %simSesultCell(w, :) = {centralized_TxPower_allWBSs_allRuns,random_TxPower_allWBSs_allRuns, centralized_CellThrought_allWBSs_allRuns, random_CellThrought_allWBSs_allRuns};

            
end

simSesult = cell2mat(simSesultCell);

ECCMultipleChannelPlots(0, n, simSesult, runtimes, xstick);


    end      

    end





%% plot averageDataOverNumOfChannels and averageDataOverNumOfChannels

    toc;

    
    
    
    
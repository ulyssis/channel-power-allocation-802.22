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


close all;
echo off;
clc;
addpath("/Users/max/Documents/git_li/channel-power-allocation-802.22/cvx");
addpath("/Library/gurobi801/mac64/matlab");
gurobi_setup;
savepath;

runtimes =  10 ;  % number of simulation run
    n = 10;    % number of WBS
    c = 3;     % number of channels, remeber to modify cvx_statusMsg whose length should be c. 
    m = c;     % number of primary users, with the same number of channels 
    delta = 1*10.^(-12);   % Noise;untitled.eps
    lengthSide = 60000;
    infBound = 1*10.^(-7);     % The interfernce threshold on PU contour  
    TVpower = 0;
    pathlossfactor = 2;    
    miniP = 4; % 36dbm, the minmum power for users
    maxP = 40; % 46dbm
    nET = 10;  % number of endterminals in each WBS 
    s = 8; % set standard deviation
    coverage = lengthSide/4/2 * 0.7; % the maximal distance away from the WBS, whihc a terminal can have 
                                     % This value should consider SUcellRadius.
    eta= 1; % the discount of the sum of interference from different WBSs, to represent the interference on the measurement point 

    tic;

    PMiu = 1/100;
    POpertation = 1/PMiu;
    
    SUcellRadius = 3000;
    schemeIIEnabled = 1;
    
%for SUcellRadius = 0:1000:1000 % 1000:1000:7000
        
        
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


for run = 1: runtimes % the number of simulations
    
%         %% make sure that both LP and CVX are feasible
%         linprogWork = -2;
%         max_cvx_statusMsg =2;
%         while (linprogWork <= 0 && max_cvx_statusMsg > 1)
%             [posSU, posET, posTVContor, Gtilde, GtildeETsSUs, GtildeAll] = geoinfo(n, m, nET, lengthSide, coverage, SUcellRadius, pathlossfactor, s);            
%             [P_LP, linprogWork] = maximalPowerPlanningLP(n, m, infBound, GtildeAll, miniP, maxP);      % The fixed power levels on all channels for every node.
%             [P_CVX, max_cvx_statusMsg] = maximalPowerPlanningCVX(n, m, infBound, GtildeAll, miniP, maxP);
%         end
%         
%         lp = sum(P_LP,2);
%         lp_container = [lp_container, lp];


    %%---------- LP ---------------    
    %   decide the maximal transmission power by solving the Linear problem with matlab.
    %   'while' loop is used to generate solution feasible topologies  
    
%         linprogWork = -2;
%         while (linprogWork <= 0)
%             [posSU, posET, posTVContor, Gtilde, GtildeETsSUs, GtildeAll] = geoinfo(n, m, nET, lengthSide, coverage, SUcellRadius, pathlossfactor, s);            
%             [P_LP, linprogWork] = maximalPowerPlanningLP(n, m, infBound, GtildeAll, miniP, maxP);      % The fixed power levels on all channels for every node.
%         end
%             lp = sum(P_LP,2);
%             lp_container = [lp_container, lp];

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


% %% plot the distibution of stations and nodes
%     [posET] = ETlocation(n, nET, lengthSide, coverage);
%     plotlocation(n, m, lengthSide, posSU, posET, posTVContor);
%     %plotMaximalPower(P, n, c);


% %% run channel assignment scheme I and comparison schemes
%         [utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
%             SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, ...
%             B_random, B_cat, B_case, B_optimization, B_noregret, B_PotentialGame, RetGUROBI_FCC] ...
%             = runSchemes(run, P_CVX, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory, powerHistory, ...
%             averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
%             SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, ...
%                 schemeIIEnabled, PMiu, POpertation, infBound, RetGUROBI_FCC);

%% run channel assignment scheme II

        [RetGUROBI_FCC] = runSchemes(run, P_CVX, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory, powerHistory, ...
            averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
            SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_optimization_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_optimization_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, ...
                schemeIIEnabled, PMiu, POpertation, infBound, RetGUROBI_FCC);


end

plotLog = SUcellRadius/100+10;

%% plot schemes

        printplotsCAschemes2power(plotLog, n, nET, powerHistory, averageSinrHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_optimization_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container);        
%         printplotsWorst20_CAschemes2power(plotLog, n, worstSINR_random_container,worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, worstSINR_optimization_container);
        
%         % plot figures for only maximal power decision
%         worst20(n, worstSINR_random_container,worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, worstSINR_optimization_container);
%         worst20(n, worstSINR_random_container2,worstSINR_cat_container2, worstSINR_case_container2, worstSINR_noregret_container2, worstSINR_optimization_container2);
%end




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

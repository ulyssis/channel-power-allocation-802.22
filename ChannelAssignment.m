%---------------------------------------------------------
%--     Channel seclection in cognitive radio networks.
%---------------------------------------------------------
% distributed greedy search, feedback infomation is used in the utility
% each user looks for the channel to minimize its utility.
% 4 schemes: 'random', 'update' and 'selfishUpdate', 'noregret'. 
% Di Li & James Gross, "DiCAPS" step 1, Written by Di Li 28/10/11

% utility fucntion: U = f_i/P_i + \sum (f_{ij}/P_j), where f is the
% interferece, P_i is the transmission power, Channel c is ommited here for clearity.

% run four distributed schemes under different ways of deciding the maximal
% transmission power!

% to use cvx:
% cd /home/li/work/tools/cvx/
% cvx_setup

% to use gurobi
% cd /home/li/work/gurobi/gurobi602/linux64/matlab
% gurobi_setup

% if lindo is called, the following step is needed to automize the running.
% cd /home/li/work/dev/DiCAPS/2011OtcGameOptReform/codes_cameraready
% mex readMatrixB.cpp 
% mex readMatrixB_jointPowerChannelAllocation.cpp

% clear all;
close all;
echo off;
clc;

runtimes =  50 ;  % number of simulation run
    n = 16;    % number of WBS
    c = 4;     % number of channels, remeber to modify cvx_statusMsg whose length should be c. 
    m = c;     % number of primary users, with the same number of channels 
    delta = 1*10.^(-12);   % Noise;untitled.eps
    lengthSide = 60000;
    SUcellRadius = 7000; % the radius of WRAN cell 3
    infBound = 1*10.^(-7);     % The interfernce threshold on PU contour    
%     infBound = 2*10.^(-7);     % The interfernce threshold on PU contour
    TVpower = 0;    
    pathlossfactor = 2;    
    miniP = 4; % 36dbm, the minmum power for users
    maxP = 40; % 46dbm
    nET = 50;   % number of endterminals in each WBS 
    s = 8; % set standard deviation
    coverage = lengthSide/4/2; % the radius of a WBS cell
    eta= 1; % the discount ofrapist the sum 5.7580of interference from different WBSs, to represent the interference on the measurement point 

    tic;

    
    
    
for SUcellRadius = 1000:1000:7000
        
        
utilityHistory=[];
powerHistory=[];      % record thress items in 10 runs of simulation
averageSinrHistory = [];
averageStdHistory = [];

SINR_ETs_random_container = [];     
SINR_ETs_whitecat_container = [];
SINR_ETs_whitecase_container = [];    
SINR_ETs_noregret_container = [];  
SINR_ETs_PotentialGame_container = [];  
SINR_ETs_lindo_container = [];
SINR_ETs_lindo2_container = [];

fair_random_container = [];     
fair_cat_container = [];     
fair_case_container = [];     
fair_lindo_container =[];
fair_noregret_container = [];  
fair_PotentialGame_container = [];  

worstSINR_random_container = [];
worstSINR_cat_container = [];     
worstSINR_case_container = [];     
worstSINR_lindo_container = [];     
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

%   ++++
utilityHistory2=[];
powerHistory2 = [];      % record thress items in 10 runs of simulation
averageSinrHistory2 = [];
averageStdHistory2 = [];

SINR_ETs_random_container2 = [];     
SINR_ETs_whitecat_container2 = [];
SINR_ETs_whitecase_container2 = [];    
SINR_ETs_noregret_container2 = [];  
SINR_ETs_PotentialGame_container2 = [];  
SINR_ETs_lindo_container2 = [];
SINR_ETs_lindo2_container2 = [];

fair_random_container2 = [];     
fair_cat_container2 = [];     
fair_case_container2 = [];    
fair_lindo_container2 =[];
fair_noregret_container2 = [];  
fair_PotentialGame_container2 = [];

worstSINR_random_container2 = [];
worstSINR_cat_container2 = [];     
worstSINR_case_container2 = [];     
worstSINR_noregret_container2 = [];   
worstSINR_lindo_container2 = [];
worstSINR_PotentialGame_container2 = [];  

convergenceStepWhitecat2 = zeros(1, runtimes);
convergenceStepWhitecase2 = zeros(1, runtimes);
convergenceStepNoregret2 = zeros(1, runtimes);
convergenceStepPotentialGame2 = zeros(1, runtimes);

SINRvarianceWhitecat_container2 = zeros(1, runtimes);
SINRvarianceWhitecase_container2 = zeros(1, runtimes);
SINRvarianceNoregret_container2 = zeros(1, runtimes);
SINRvariancePotentialGame_container2 = zeros(1, runtimes);

B_random =[];
B_cat =[];
B_case=[];
B_lindo=[];
B_noregret=[];

lp_container=[];
cvx_container=[];


% after power allocation
B_powerAllocation =[];
utilityHistory_pa=[];
powerHistory_pa=[];
averageSinrHistory_pa=[];
averageStdHistory_pa=[];
SINR_ETs_random_container_pa=[];
SINR_ETs_whitecat_container_pa=[];
SINR_ETs_whitecase_container_pa=[];
SINR_ETs_lindo_container_pa=[];
SINR_ETs_noregret_container_pa=[];
SINR_ETs_PotentialGame_container_pa=[];
SINR_ETs_LindoCAPA_container_pa=[];









for run = 1: runtimes % the number of simulations
        %   % make sure that both LP and CVX are feasible
        linprogWork = -2;
        max_cvx_statusMsg =2;
        while (linprogWork <= 0 && max_cvx_statusMsg > 1)
            [posSU, posET, posTVContor, Gtilde, GtildeETsSUs, GtildeAll] = geoinfo(n, m, nET, lengthSide, coverage, SUcellRadius, pathlossfactor, s);            
            [P_LP, linprogWork] = maximalPowerPlanningLP(n, m, infBound, GtildeAll, miniP, maxP);      % The fixed power levels on all channels for every node.
            [P_CVX, max_cvx_statusMsg] = maximalPowerPlanningCVX(n, m, infBound, GtildeAll, miniP, maxP);
        end
        
        lp = sum(P_LP,2);
        lp_container = [lp_container, lp];

        
    %  %---------- LP ---------------    
    %   decide the maximal transmission power by solving the Linear problem with matlab.
    %   'while' loop is used to generate solution feasible topologies  
    
        linprogWork = -2;
        while (linprogWork <= 0)
            [posSU, posET, posTVContor, Gtilde, GtildeETsSUs, GtildeAll] = geoinfo(n, m, nET, lengthSide, coverage, SUcellRadius, pathlossfactor, s);            
            [P_LP, linprogWork] = maximalPowerPlanningLP(n, m, infBound, GtildeAll, miniP, maxP);      % The fixed power levels on all channels for every node.
        end
            lp = sum(P_LP,2);
            lp_container = [lp_container, lp];

    %   %---------- cvx ---------------
    %    decide the maximal transmission power by solving the convex problem with cvx 
    %    while loop is used to generate solution feasible topologies
        max_cvx_statusMsg =2;
        while(max_cvx_statusMsg > 1) % contains cvx_optval which is the value of the objective function 
            [posSU, posET, posTVContor, Gtilde, GtildeETsSUs, GtildeAll] = geoinfo(n, m, nET, lengthSide, coverage, SUcellRadius, pathlossfactor, s);            
            [P_CVX, max_cvx_statusMsg] = maximalPowerPlanningCVX(n, m, infBound, GtildeAll, miniP, maxP);
        end
            cvx=sum(P_CVX,2);
            cvx_container = [cvx_container, cvx];


%   % plot the distibution of stations and nodes
%     [posET] = ETlocation(n, nET, lengthSide, coverage);
%     plotlocation(n, m, lengthSide, posSU, posET, posTVContor);
%     plotMaximalPower(P, n, c);

%%%% run channel assignment schemes
        [utilityHistory, powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
            SINR_ETs_lindo_container, SINR_ETs_lindo2_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_lindo_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_lindo_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container, ...
            B_random, B_cat, B_case, B_lindo, B_noregret, B_PotentialGame, B_lindoCAPA] ...%snrRatio_random, snrRatio_dica, snrRatio_self, snrRatio_noregret]... % to help function printplots4schemes
            = runSchemes(run, P_CVX, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory, powerHistory, ...
            averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, ...
            SINR_ETs_lindo_container, SINR_ETs_lindo2_container, SINR_ETs_noregret_container, SINR_ETs_PotentialGame_container, fair_random_container, fair_cat_container, fair_case_container, fair_lindo_container, fair_noregret_container, fair_PotentialGame_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_lindo_container, worstSINR_noregret_container, worstSINR_PotentialGame_container, convergenceStepWhitecat, convergenceStepWhitecase, convergenceStepNoregret, convergenceStepPotentialGame, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container, SINRvariancePotentialGame_container);
        
        
% %%% run power allocation after channel allocation is completed.
% %%% input: channel allocation. 
% %%% system parameters
%     scheme = [{'random'}, {'whitecat'}, {'whitecase'}, {'lindo'}, {'noregret'}];
%     for schemeIndex = 1: 5;
% 
%         if (strcmpi(scheme(schemeIndex), 'random'))
%             B = B_random;
%         end
% 
%         if (strcmpi(scheme(schemeIndex), 'whitecat'))
%             B = B_cat;
%         end
% 
%         if (strcmpi(scheme(schemeIndex), 'whitecase'))
%             B = B_case;
%         end
% 
%         if (strcmpi(scheme(schemeIndex), 'lindo'))
%             B = B_lindo;
%         end
% 
%         if (strcmpi(scheme(schemeIndex), 'noregret'))
%             B = B_noregret;
%         end
% 
%         %  do power allocation after channel allocation
%         [B_powerAllocation] = powerAllocation(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
%         [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B_powerAllocation, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
%         [SINR_ETs_pa, worstSINR_pa, fair_pa] = SINR_ETs_cellReSelection(B_powerAllocation, n, GtildeETsSUs, nET, TVpower, delta);
% 
%         if (strcmpi(scheme(schemeIndex), 'random'))
%             random_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
%             SINR_ETs_random_container_pa = [SINR_ETs_random_container_pa, SINR_ETs_pa];
%         end
% 
%         if (strcmpi(scheme(schemeIndex), 'whitecat'))
%             whitecat_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
%             SINR_ETs_whitecat_container_pa = [SINR_ETs_whitecat_container_pa, SINR_ETs_pa];            
%         end
% 
%         if (strcmpi(scheme(schemeIndex), 'whitecase'))
%             whitecase_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
%             SINR_ETs_whitecase_container_pa = [SINR_ETs_whitecase_container_pa, SINR_ETs_pa];            
%             
%         end
% 
%         if (strcmpi(scheme(schemeIndex), 'lindo'))
%             lindo_Perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
%             SINR_ETs_lindo_container_pa = [SINR_ETs_lindo_container_pa, SINR_ETs_pa];            
%         end
% 
%         if (strcmpi(scheme(schemeIndex), 'noregret'))
%             noregret_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
%             SINR_ETs_noregret_container_pa = [SINR_ETs_noregret_container_pa, SINR_ETs_pa];           
%         end
%         
%     end
%     
%     % two joint-power/channel allocation
%         scheme = [{'PotentialGame'}, {'Lindo-CA+PA'}];
%     for schemeIndex = 1: 2;
%         if (strcmpi(scheme(schemeIndex), 'PotentialGame'))
%             B = B_PotentialGame;
%         end
%         
%         if (strcmpi(scheme(schemeIndex), 'Lindo-CA+PA'))
%             B = B_lindoCAPA;
%         end
%         
%         [sumUtility, averageI, averageP, averageSINR, stdSINR] = obtainPerformance(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor);
%         [SINR_ETs_pa, worstSINR_pa, fair_pa] = SINR_ETs_cellReSelection(B_powerAllocation, n, GtildeETsSUs, nET, TVpower, delta);
% 
%         
%         if (strcmpi(scheme(schemeIndex), 'PotentialGame'))
%             PotentialGame_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
%             SINR_ETs_PotentialGame_container_pa = [SINR_ETs_PotentialGame_container_pa, SINR_ETs_pa];           
%         end        
%         
%         if (strcmpi(scheme(schemeIndex), 'Lindo-CA+PA'))
%             LindoCAPA_perf = [sumUtility, averageI, averageP, averageSINR, stdSINR];
%             SINR_ETs_LindoCAPA_container_pa = [SINR_ETs_LindoCAPA_container_pa, SINR_ETs_pa];           
%         end            
%         
%         
%     end
%     
%     
%     
%         utilityHistory_pa(:, run) = [random_perf(1); lindo_Perf(1); whitecat_perf(1); whitecase_perf(1); noregret_perf(1); PotentialGame_perf(1) ; LindoCAPA_perf(1)];
%         powerHistory_pa(:, run) = [random_perf(3); lindo_Perf(3); whitecat_perf(3); whitecase_perf(3); noregret_perf(3); PotentialGame_perf(3); LindoCAPA_perf(3)];
%         averageSinrHistory_pa(:, run) = [random_perf(4); lindo_Perf(4); whitecat_perf(4); whitecase_perf(4); noregret_perf(4); PotentialGame_perf(4); LindoCAPA_perf(4)];
%         averageStdHistory_pa(:, run) = [random_perf(5); lindo_Perf(5); whitecat_perf(5); whitecase_perf(5); noregret_perf(5); PotentialGame_perf(5); LindoCAPA_perf(5)];  

%% plot figures for the performances after power allocation
[utilityHistory2, powerHistory2, averageSinrHistory2, averageStdHistory2, SINR_ETs_random_container2, SINR_ETs_whitecat_container2, SINR_ETs_whitecase_container2, ...
            SINR_ETs_lindo_container2, SINR_ETs_lindo2_container2, SINR_ETs_noregret_container2, SINR_ETs_PotentialGame_container2, fair_random_container2, fair_cat_container2, fair_case_container2, fair_lindo_container2, fair_noregret_container2, fair_PotentialGame_container2, worstSINR_random_container2, worstSINR_cat_container2, worstSINR_case_container2, worstSINR_lindo_container2, worstSINR_noregret_container2, worstSINR_PotentialGame_container2, convergenceStepWhitecat2, convergenceStepWhitecase2, convergenceStepNoregret2, convergenceStepPotentialGame2, SINRvarianceWhitecat_container2, SINRvarianceWhitecase_container2, SINRvarianceNoregret_container2, SINRvariancePotentialGame_container2, ...
            B_random2, B_cat2, B_case2, B_lindo2, B_noregret2, B_PotentialGame2, B_lindoCAPA2] ...%snrRatio_random, snrRatio_dica, snrRatio_self, snrRatio_noregret]... % to help function printplots4schemes
            = runSchemes(run, P_CVX, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, utilityHistory2, powerHistory2, ...
            averageSinrHistory2, averageStdHistory2, SINR_ETs_random_container2, SINR_ETs_whitecat_container2, SINR_ETs_whitecase_container2, ...
            SINR_ETs_lindo_container2, SINR_ETs_lindo2_container2, SINR_ETs_noregret_container2, SINR_ETs_PotentialGame_container2, fair_random_container2, fair_cat_container2, fair_case_container2, fair_lindo_container2, fair_noregret_container2, fair_PotentialGame_container2, worstSINR_random_container2, worstSINR_cat_container2, worstSINR_case_container2, worstSINR_lindo_container2, worstSINR_noregret_container2, worstSINR_PotentialGame_container2, convergenceStepWhitecat2, convergenceStepWhitecase2, convergenceStepNoregret2, convergenceStepPotentialGame2, SINRvarianceWhitecat_container2, SINRvarianceWhitecase_container2, SINRvarianceNoregret_container2, SINRvariancePotentialGame_container);
end

plotLog = SUcellRadius/100;
% %%%%%----channel allocation schemes, two power formulations   ------------
%         % plots of the 4 schemes under two kinds of power decision.
        printplotsCAschemes2power(plotLog, n, powerHistory, averageSinrHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_noregret_container, SINR_ETs_lindo_container, powerHistory2, averageSinrHistory2, SINR_ETs_random_container2, SINR_ETs_whitecat_container2, SINR_ETs_whitecase_container2, SINR_ETs_noregret_container2, SINR_ETs_lindo_container2);        
        printplotsWorst20_CAschemes2power(plotLog, n, worstSINR_random_container,worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, worstSINR_lindo_container, worstSINR_random_container2,worstSINR_cat_container2, worstSINR_case_container2, worstSINR_noregret_container2, worstSINR_lindo_container2);
        
%         % plot figures for only maximal power decision
%         worst20(n, worstSINR_random_container,worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, worstSINR_lindo_container);
%         worst20(n, worstSINR_random_container2,worstSINR_cat_container2, worstSINR_case_container2, worstSINR_noregret_container2, worstSINR_lindo_container2);



end

%%%-------This function generates obsolete plot--------
% %       the runSchemes function should add "snrRatio_random, snrRatio_dica,
% %       snrRatio_self, snrRatio_noregret" to the end of its output
% %       parameters
%          printplots4schemes(powerHistory, averageSinrHistory, averageStdHistory, snrRatio_random, snrRatio_dica, snrRatio_self, snrRatio_noregret);
%%%-------This function generates obsolete plot--------



% %%%-------plots: joint channel and power allocation -----------
% %%% Running 5 schemes requires results form Lindo, we call the c++ program
% %%% from matlab:
% %%%           mex –setup   % choose complier
% %%%           mex readMatrixB.cpp
% %%%           mex readMatrixB_jointPowerChannelAllocation.cpp
% %%%-------power consumption, qusai SINR, and sinr on terminals.
% plotTag = 0;
%          printplots5schemes(plotTag, n, utilityHistory, powerHistory, averageSinrHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, SINR_ETs_lindo_container, SINR_ETs_noregret_container);
% 
%         %printPerformanceETs5(SINR_ETs_random, SINR_ETs_whitecat, SINR_ETs_whitecase, SINR_ETs_lindo, SINR_ETs_noregret, fair_random_container, fair_cat_container, fair_case_container, fair_lindo, fair_noregret_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, fair_lindo, worstSINR_noregret_container);
% %         printplots5schemes(powerHistory, averageSinrHistory);
% 
% % plots after power allocation
% % 5 schemes + 2 joint schemes
% plotTag = plotTag+1;
%          printplots7schemes(plotTag, n, utilityHistory_pa, powerHistory_pa, averageSinrHistory_pa, SINR_ETs_random_container_pa, SINR_ETs_whitecat_container_pa, SINR_ETs_whitecase_container_pa, SINR_ETs_lindo_container_pa, SINR_ETs_noregret_container_pa, SINR_ETs_PotentialGame_container_pa, SINR_ETs_LindoCAPA_container_pa);








        


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
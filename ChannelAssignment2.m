% WBSs' locations are fixed, random update sequence!

% is able to run 4 distributed schemes, besides, lindo needs manual assistance

clear all;
close all;
echo off;
clc;

runtimes = 2  ;  % number of simulation run
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


powerHistory=[];      % record thress items in 10 runs of simulation
averageSinrHistory = [];
averageStdHistory = [];

SINR_ETs_random_container = [];     
SINR_ETs_whitecat_container = [];
SINR_ETs_whitecase_container = [];
sinr_ETs_lindo_container = [];
SINR_ETs_lindo2_container = [];
SINR_ETs_noregret_container = [];  

fair_random_container = [];     
fair_cat_container = [];     
fair_case_container = [];     
fair_noregret_container = [];  

worstSINR_random_container = [];
worstSINR_cat_container = [];     
worstSINR_case_container = [];     
worstSINR_noregret_container = [];     

convergenceStepWhitecat = zeros(1, runtimes);
convergenceStepWhitecase = zeros(1, runtimes);
convergenceStepNoregret = zeros(1, runtimes);

SINRvarianceWhitecat_container = zeros(1, runtimes);
SINRvarianceWhitecase_container = zeros(1, runtimes);
SINRvarianceNoregret_container = zeros(1, runtimes);

%   ++++
powerHistory2 = [];      % record thress items in 10 runs of simulation
averageSinrHistory2 = [];
averageStdHistory2 = [];

SINR_ETs_random_container2 = [];     
SINR_ETs_whitecat_container2 = [];
SINR_ETs_whitecase_container2 = [];    
SINR_ETs_noregret_container2 = [];  

fair_random_container2 = [];     
fair_cat_container2 = [];     
fair_case_container2 = [];     
fair_noregret_container2 = [];  

worstSINR_random_container2 = [];
worstSINR_cat_container2 = [];     
worstSINR_case_container2 = [];     
worstSINR_noregret_container2 = [];     

convergenceStepWhitecat2 = zeros(1, runtimes);
convergenceStepWhitecase2 = zeros(1, runtimes);
convergenceStepNoregret2 = zeros(1, runtimes);

SINRvarianceWhitecat_container2 = zeros(1, runtimes);
SINRvarianceWhitecase_container2 = zeros(1, runtimes);
SINRvarianceNoregret_container2 = zeros(1, runtimes);


lp_container=[];
cvx_container=[];



% use cvx, fixed TV contours
       max_cvx_statusMsg =2;
       while (max_cvx_statusMsg > 1)% find feasible locations to apply LP and CVX
            [posSU, posET, posTVContor, Gtilde, GtildeETsSUs, GtildeAll] = geoinfo(n, m, nET, lengthSide, coverage, SUcellRadius, pathlossfactor, s);            
            [P_CVX, max_cvx_statusMsg] = maximalPowerPlanningCVX(n, m, infBound, GtildeAll, miniP, maxP);
       end

       cvx=sum(P_CVX,2);
       cvx_container = [cvx_container, cvx];



% [posET] = ETlocation(n, nET, lengthSide, coverage);
for run = 1: runtimes % the number of simulations

%       plotlocation(n, m, lengthSide, posSU, posET, posTVContor);
%     plotMaximalPower(P, n, c);

            [powerHistory, averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, sinr_ETs_lindo_container,...
            SINR_ETs_lindo2_container, SINR_ETs_noregret_container, fair_random_container, fair_cat_container, fair_case_container, fair_noregret_container, worstSINR_random_container, ...
            worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, convergenceStepWhitecat, convergenceStepWhitecase, ...
            convergenceStepNoregret, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container]...
            = runSchemes(run, infBound, P_CVX, Gtilde, GtildeETsSUs, n, c, m, nET, GtildeAll, TVpower, SUcellRadius, delta, pathlossfactor, eta, powerHistory, ...
            averageSinrHistory, averageStdHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, sinr_ETs_lindo_container,...
            SINR_ETs_lindo2_container, SINR_ETs_noregret_container, fair_random_container, fair_cat_container, fair_case_container, fair_noregret_container, worstSINR_random_container, ...
            worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, convergenceStepWhitecat, convergenceStepWhitecase, ...
            convergenceStepNoregret, SINRvarianceWhitecat_container, SINRvarianceWhitecase_container, SINRvarianceNoregret_container);

end
        
            %printplots5schemes(n, powerHistory, averageSinrHistory, lindo_averagePower, lindo_averageSINR, lindo2_averagePower, lindo2_averageSINR, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, sinr_ETs_lindo_container, SINR_ETs_lindo2_container, SINR_ETs_noregret_container);

            printplots5schemes(n, powerHistory, averageSinrHistory, SINR_ETs_random_container, SINR_ETs_whitecat_container, SINR_ETs_whitecase_container, sinr_ETs_lindo_container, SINR_ETs_noregret_container);


    toc;

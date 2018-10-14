%  plot performance of 4 distributed schemes 
%  output1: average Power
%  output2: average QuasiSINR
%  output3: sinr on end users


function printplotsCAschemes2power(plotLog, n, nET, powerHistory, averageSinrHistory, SINR_ETs_random, SINR_ETs_whitecat, SINR_ETs_whitecase, SINR_ETs_lindo, SINR_ETs_noregret, SINR_ETs_PotentialGame)
% the latter snrRatio_random, snrRatio_whitecat, snrRatio_whitecase,
% snrRatio_noregret are sinr on each WBS in the last run, which is a
% intersection for sinr in one run

%% Average Sinr
figure(plotLog+0);
handle1 = barweb(mean(averageSinrHistory,2)', 1.96*std(averageSinrHistory,1,2)'/sqrt(n), [], [], [], [], 'Average SINR (dB)', bone, 'y', {'Random Allocation'; 'whiteCat'; 'whiteCase'; 'Noregret Learning' ; 'optimization'; 'Potential Game'}, 2, 'plot');
set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');

%% Average Transmisson Power
figure(plotLog+1);
% use barweb function 
% dbstop on error;
handle1 = barweb(mean(powerHistory,2)', 1.96*std(powerHistory,1,2)'/sqrt(n), [], [], [], [], 'Average transmission power of each WBS (W)', bone, 'y', {'Random Allocation'; 'whiteCat'; 'whiteCase'; 'Noregret Learning' ; 'optimization'; 'Potential Game'}, 2, 'plot');
set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');

%% SINR cfg
figure (plotLog+3)
    h1 = cdfplot(SINR_ETs_random);
    set(h1, 'Color','k', 'LineStyle', '--');
    hold on
    h2 = cdfplot(SINR_ETs_whitecat);
    set(h2, 'Color','g', 'LineStyle', '-');
    hold on
    h3 = cdfplot(SINR_ETs_whitecase);
    set(h3, 'Color', 'r', 'LineStyle', ':');
    hold on
    h5 = cdfplot(SINR_ETs_noregret);
    set(h5,'Color','b',  'LineStyle', ':');
    hold on
    h4 = cdfplot(SINR_ETs_lindo);
    set(h4,'Color','m','LineStyle', '--');
    hold on
    h6 = cdfplot(SINR_ETs_PotentialGame);
    set(h6,'Color','c',  'LineStyle', '-.');
    
    axis([1 100 0 1]);
    set(gca, 'XScale', 'log');
    yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
    yticklabels({'0','10%','20%','30%','40%','50%','60%', '70%', '80%', '90%', '100%'})

    legend('Random Allocation', 'whiteCat', 'whiteCase', 'Noregret Learning', 'optimization', 'Potential Game', 'Location', 'southeast');
   
    xlabel('SINR (dB)');
    set(gca,'YLabel',[]);
    magnify;
   
%% SINR average overall all end users form all runs. should not be used.
figure (plotLog+4)
averageSINR = [mean(SINR_ETs_random), mean(SINR_ETs_whitecat), mean(SINR_ETs_whitecase), mean(SINR_ETs_noregret), mean(SINR_ETs_lindo), mean(SINR_ETs_PotentialGame)];
stdSINR = [1.96*std(SINR_ETs_random,1)/sqrt(size(SINR_ETs_random, 2)), 1.96*std(SINR_ETs_whitecat,1)/sqrt(size(SINR_ETs_whitecat, 2)), 1.96*std(SINR_ETs_whitecase,1)/sqrt(size(SINR_ETs_whitecase, 2)), 1.96*std(SINR_ETs_noregret,1)/sqrt(size(SINR_ETs_noregret, 2)), 1.96*std(SINR_ETs_lindo,1)/sqrt(size(SINR_ETs_lindo, 2)), 1.96*std(SINR_ETs_PotentialGame,1)/sqrt(size(SINR_ETs_PotentialGame, 2))];
handle1 = barweb(averageSINR, stdSINR, [], [], [], [], 'Average SINR on End Users over All Runs (dB) (not useful)', bone, 'y', {'Random Allocation'; 'whiteCat'; 'whiteCase'; 'Noregret Learning' ; 'optimization'; 'Potential Game'}, 2, 'plot');
set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');


%% SINR average on end users.

figure (plotLog+5)

nRuns = size(SINR_ETs_random, 2)/n/nET;

average_SINR_ETs_random = zeros(1, nRuns);
% CI_SINR_ETs_random = zeros(1, nRuns);
for i = 1:nRuns
    average_SINR_ETs_random(i) = mean(SINR_ETs_random((i-1)*nET*n+1: i*nET*n ));
%     CI_SINR_ETs_random(i) = 1.96*std(SINR_ETs_random((i-1)*nET*n+1: i*nET*n ),1)/sqrt(nET*n);
end

average_SINR_ETs_whitecat = zeros(1, nRuns);
for i = 1:nRuns
    average_SINR_ETs_whitecat(i) = mean(SINR_ETs_whitecat((i-1)*nET*n+1: i*nET*n ));
end

average_SINR_ETs_whitecase = zeros(1, nRuns);
for i = 1:nRuns
    average_SINR_ETs_whitecase(i) = mean(SINR_ETs_whitecase((i-1)*nET*n+1: i*nET*n ));
end
    
average_SINR_ETs_noregret = zeros(1, nRuns);
for i = 1:nRuns
    average_SINR_ETs_noregret(i) = mean(SINR_ETs_noregret((i-1)*nET*n+1: i*nET*n ));
end

average_SINR_ETs_lindo = zeros(1, nRuns);
for i = 1:nRuns
    average_SINR_ETs_lindo(i) = mean(SINR_ETs_lindo((i-1)*nET*n+1: i*nET*n ));
end

average_SINR_ETs_PotentialGame = zeros(1, nRuns);
for i = 1:nRuns
    average_SINR_ETs_PotentialGame(i) = mean(SINR_ETs_PotentialGame((i-1)*nET*n+1: i*nET*n ));
end

average_SINR_ETs = [average_SINR_ETs_random ;
    average_SINR_ETs_whitecat;
    average_SINR_ETs_whitecase;
    average_SINR_ETs_noregret;
    average_SINR_ETs_lindo;
    average_SINR_ETs_PotentialGame];


handle1 = barweb(mean(average_SINR_ETs,2)', 1.96*std(average_SINR_ETs,1,2)'/sqrt(n), [], [], [], [], 'Average SINR on End Terminals (dB)', bone, 'y', {'Random Allocation'; 'whiteCat'; 'whiteCase'; 'Noregret Learning' ; 'optimization'; 'Potential Game'}, 2, 'plot');
set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');




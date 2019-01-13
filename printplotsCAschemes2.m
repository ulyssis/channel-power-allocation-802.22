%  plot performance of Scheme II
%  output1: average Power
%  output2: average QuasiSINR
%  output3: sinr on end users


function printplotsCAschemes2(plotLog, n, nET, utilityHistoryFCC, powerHistoryFCC, averageSinrHistoryFCC, SINR_ETs_centralized_FCC_container, SINR_ETs_distributed_FCC_container)

% the latter snrRatio_random, snrRatio_whitecat, snrRatio_whitecase,
% snrRatio_noregret are sinr on each WBS in the last run, which is a
% intersection for sinr in one run

%% Average Sinr
figure(plotLog+0);
handle1 = barweb(mean(averageSinrHistoryFCC,2)', 1.96*std(averageSinrHistoryFCC,1,2)'/sqrt(n), [], [], [], [], 'Average SINR (dB)', bone, 'y', {'Scheme II centralized'; 'Scheme II distributed'}, 2, 'plot');
set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');

%% Average Transmisson Power
figure(plotLog+1);
% use barweb function 
% dbstop on error;
handle1 = barweb(mean(powerHistoryFCC,2)', 1.96*std(powerHistoryFCC,1,2)'/sqrt(n), [], [], [], [], 'Average transmission power of each WBS (W)', bone, 'y', {'Scheme II centralized'; 'Scheme II distributed'}, 2, 'plot');
set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');

%% sumUtility
figure(plotLog+2);

% delete the all zero colume.
utilityHistoryFCC = utilityHistoryFCC(:, any(utilityHistoryFCC,1));

handle1 = barweb(mean(utilityHistoryFCC,2)', 1.96*std(utilityHistoryFCC,1,2)'/sqrt(n), [], [], [], [], 'Sum of Utility of All WBSs', bone, 'y', {[' ']}, 2, 'plot'); %['Scheme II centralized'; 'Scheme II distributed']
set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');

figure(plotLog*10+2);

% delete the all zero colume.
utilityHistoryFCC = utilityHistoryFCC(:, any(utilityHistoryFCC,1));

%x = ['Scheme II centralized', 'Scheme II distributed'];
data = mean(utilityHistoryFCC,2)';
errhigh = 1.96*std(utilityHistoryFCC,1,2)'/sqrt(n);
errlow  = 1.96*std(utilityHistoryFCC,1,2)'/sqrt(n);

bar(x,data);  
xticks([1 2])
xticklabels({'Scheme II centralized', 'Scheme II distributed'})
ylabel('Sum of Utility of All WBSs');

hold on

er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

set(gca,'fontsize',10);
hold off

%% SINR cfg
figure (plotLog+3)
    h1 = cdfplot(SINR_ETs_centralized_FCC_container);
    set(h1, 'Color','k', 'LineStyle', '--');
    hold on
    h2 = cdfplot(SINR_ETs_distributed_FCC_container);
    set(h2, 'Color','g', 'LineStyle', '-');

    
    axis([1 100 0 1]);
    set(gca, 'XScale', 'log');
    yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
    yticklabels({'0','10%','20%','30%','40%','50%','60%', '70%', '80%', '90%', '100%'})

    legend('Scheme II centralized', 'Scheme II distributed', 'Location', 'southeast');
   
    xlabel('SINR (dB)');
    set(gca,'YLabel',[]);
    magnify;
   
% %% SINR average overall all end users form all runs. should not be used.
% figure (plotLog+4)
% averageSINR = [mean(SINR_ETs_centralized_FCC_container), mean(SINR_ETs_distributed_FCC_container)];
% stdSINR = [1.96*std(SINR_ETs_centralized_FCC_container,1)/sqrt(size(SINR_ETs_centralized_FCC_container, 2)), 1.96*std(SINR_ETs_distributed_FCC_container,1)/sqrt(size(SINR_ETs_distributed_FCC_container, 2))];
% handle1 = barweb(averageSINR, stdSINR, [], [], [], [], 'Average SINR on End Users over All Runs (dB) (not useful)', bone, 'y', {'Scheme II centralized'; 'Scheme II distributed'}, 2, 'plot');
% set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');


%% SINR average on end users.

figure (plotLog+5)

nRuns = size(SINR_ETs_centralized_FCC_container, 2)/n/nET;

average_SINR_ETs_centralized_FCC_container = zeros(1, nRuns);
% CI_SINR_ETs_centralized_FCC_container = zeros(1, nRuns);
for i = 1:nRuns
    average_SINR_ETs_centralized_FCC_container(i) = mean(SINR_ETs_centralized_FCC_container((i-1)*nET*n+1: i*nET*n ));
%     CI_SINR_ETs_centralized_FCC_container(i) = 1.96*std(SINR_ETs_centralized_FCC_container((i-1)*nET*n+1: i*nET*n ),1)/sqrt(nET*n);
end

average_SINR_ETs_distributed_FCC_container = zeros(1, nRuns);
for i = 1:nRuns
    average_SINR_ETs_distributed_FCC_container(i) = mean(SINR_ETs_distributed_FCC_container((i-1)*nET*n+1: i*nET*n ));
end



average_SINR_ETs = [average_SINR_ETs_centralized_FCC_container ;
    average_SINR_ETs_distributed_FCC_container];


handle1 = barweb(mean(average_SINR_ETs,2)', 1.96*std(average_SINR_ETs,1,2)'/sqrt(n), [], [], [], [], 'Average SINR on End Terminals (dB)', bone, 'y', {'Scheme II centralized'; 'Scheme II distributed'}, 2, 'plot');
set(handle1.legend,'Location','southeast', 'FontSize', 10, 'Color', 'R');




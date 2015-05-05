%function    printPerformanceETs(SINR_ETs_random, SINR_ETs_whitecat, SINR_ETs_whitecase, SINR_ETs_noregret)
% figure(5), figure(6), figure(7)
function    printPerformanceETs4(SINR_ETs_random, SINR_ETs_whitecat, SINR_ETs_whitecase, SINR_ETs_noregret, fair_random_container, fair_cat_container, fair_case_container, fair_noregret_container, worstSINR_random_container, worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container)

hold off;
figure (5)
% h = cdfplot(SINR_ETs_random; SINR_ETs_whitecat; SINR_ETs_whitecase; SINR_ETs_noregret);

    h1 = cdfplot(SINR_ETs_random);
    set(h1, 'Color','k', 'LineStyle', '-');
    hold on
    h2 = cdfplot(SINR_ETs_whitecat);
    set(h2, 'Color','g', 'LineStyle', '--');
    hold on
    h3 = cdfplot(SINR_ETs_whitecase);
    set(h3, 'Color', 'r', 'LineStyle', ':');
    hold on
%     h4 = cdfplot(SINR_ETs_lindo);420
%     set(h4,'Color','m');
%     hold on
    h5 = cdfplot(SINR_ETs_noregret);
    set(h5,'Color','b',  'LineStyle', '-.');


    % Convert y-axis values to percentage values by multiplication
    a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
    % Create a vector of '%' signs
    pct = char(ones(size(a,1),1)*'%'); 
    % Append the '%' signs after the percentage values
    new_yticks = [char(a),pct];
    % 'Reflect the changes on the plot
    set(gca,'yticklabel',new_yticks);

    
    set (gcf,'Position',[232 246 560 300]); % size of plot
    set(gca, 'Position', [.13 .17 .80 .74]); % the percentage of X, Y axis in the plot
    set(gca, 'XScale', 'log');
    legend('random', 'whitecat','whitecase', 'noregret', 'Location', 'southeast');
    title('Cumulative distribution of SINR on end terminals for different schemes');
    axis([1 100 0 1]);
%     hold on
    h5 = cdfplot(SINR_ETs_noregret);
    set(h5,'Color','b',  'LineStyle', '-.');


    % Convert y-axis values to percentage values by multiplication
    a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
    % Create a vector of '%' signs
    pct = char(ones(size(a,1),1)*'%'); 
    % Append the '%' signs after the percentage values
    new_yticks = [char(a),pct];
    % 'Reflect the changes on the plot
    set(gca,'yticklabel',new_yticks);

    
    set (gcf,'Position',[232 246 560 300]); % size of plot
    set(gca, 'Position', [.13 .17 .80 .74]); % the percentage of X, Y axis in the plot
    set(gca, 'XScale', 'log');
    legend('random', 'whitecat', 'whitecase', 'noregret', 'Location', 'southeast');
    title('Cumulative distribution of SINR on end terminals for different schemes');
    axis([1 100 0 1]);
    xlabel('SINR (db)');
%     set(gca, 'XTick', 1:4, 'XTickLabel', labels);            
%     ylabel('%')
%     s10*log10(et(gca, 'YTickLabel', percentsz);    
    


figure (6)
% Average value of the 20% worst SINR

random_mean_worst = mean(worstSINR_random_container);

whitecat_mean_worst = mean(worstSINR_cat_container);

whitecase_mean_worst = mean(worstSINR_case_container);

noregret_mean_worst = mean(worstSINR_noregret_container);

labels = {'rand' 'whitecat' 'whitecase'  'noreg'};
bar((1:4), [random_mean_worst, whitecat_mean_worst, whitecase_mean_worst, noregret_mean_worst]);
title('Average value of the 20% worst SINR');
set(gca, 'XTick', 1:4, 'XTickLabel', labels);



figure(7) 
% standard divation on the number of ETs in different cells
labels = {'rand' 'whitecat' 'whitecase'  'noreg'};
bar((1:4), [mean(fair_random_container), mean(fair_cat_container), mean(fair_case_container), mean(fair_noregret_container)]);
title(['standard deviation of the cell size',sprintf('\n'),' (the number of served end users) of all WBSs']);
set(gca, 'XTick', 1:4, 'XTickLabel', labels);

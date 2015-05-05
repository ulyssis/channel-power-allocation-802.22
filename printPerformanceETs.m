%function    printPerformanceETs(SINR_ETs_random, SINR_ETs_whitecat, SINR_ETs_whitecase, SINR_ETs_noregret)
function    printPerformanceETs(SINR_ETs_random, SINR_ETs_whitecat, SINR_ETs_whitecase, SINR_ETs_lindo, SINR_ETs_noregret)

hold off;
figure (5)
% h = cdfplot(SINR_ETs_random; SINR_ETs_whitecat; SINR_ETs_whitecase; SINR_ETs_noregret);

    h1 = cdfplot(SINR_ETs_random);
    set(h1,'Color','k');
    hold on
    h2 = cdfplot(SINR_ETs_whitecat);
    set(h2,'Color','g');
    hold on
    h3 = cdfplot(SINR_ETs_whitecase);
    set(h3,'Color','r');
    hold on
    h4 = cdfplot(SINR_ETs_lindo);
    set(h4,'Color','m');
    hold on
    h5 = cdfplot(SINR_ETs_noregret);
    set(h5,'Color','b');
    set(gca, 'XScale', 'log');
    legend('random', 'whitecat','whitecase', 'Lindo', 'noregret', 'Location', 'southeast');
    title('Cumulative distribution of SINR on end terminals')
    xlabel('SINR (db)')
%     set(gca, 'XTick', 1:4, 'XTickLabel', labels);            
    ylabel('%')
%     set(gca, 'YTickLabel', percentsz);  


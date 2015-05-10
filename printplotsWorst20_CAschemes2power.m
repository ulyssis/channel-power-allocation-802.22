function  printplotsWorst20_CAschemes2power(plotLog, n, worstSINR_random_container,worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, worstSINR_lindo_container, worstSINR_random_container2,worstSINR_cat_container2, worstSINR_case_container2, worstSINR_noregret_container2, worstSINR_lindo_container2)

    
  worstSINR = [worstSINR_random_container;worstSINR_cat_container; worstSINR_case_container; worstSINR_noregret_container; worstSINR_lindo_container];
  worstSINR2 = [worstSINR_random_container2;worstSINR_cat_container2; worstSINR_case_container2; worstSINR_noregret_container2; worstSINR_lindo_container2];
  
  
  figure (plotLog+5)
handle1=barweb([mean(worstSINR2,2)';mean(worstSINR,2)'], [1.96*std(worstSINR2,1,2)'/sqrt(n); 1.96*std(worstSINR,1,2)'/sqrt(n)], [], {'linear optimization formulation'; 'convex optimization formulation'}, [], [], 'worst 20\% SINR on end users (dB)', bone, 'y', {'random allocation'; 'whiteCat'; 'whiteCase'; 'noregret learning' ; 'optimization'}, 2, 'plot')
  set(handle1.legend,'Location','north', 'FontSize', 10, 'Color', 'R');
    set(handle1.legend, 'xcolor', 'w','ycolor', 'w');
% output: worst 20%

function worst20(n, worstSINR_random_container,worstSINR_cat_container, worstSINR_case_container, worstSINR_noregret_container, worstSINR_lindo_container)


figure (6) % The average value of the worst 20% SINR

barweb([mean(worstSINR_random_container) mean(worstSINR_cat_container) mean(worstSINR_case_container) mean(worstSINR_noregret_container) mean(worstSINR_lindo_container)], [1.96*std(worstSINR_random_container,1)/sqrt(n) 1.96*std(worstSINR_cat_container,1)/sqrt(n) 1.96*std(worstSINR_case_container,1)/sqrt(n) 1.96*std(worstSINR_noregret_container,1)/sqrt(n)  1.96*std(worstSINR_lindo_container,1)/sqrt(n)], [], [], [], [], 'SINR(dB)', bone, 'y', {'random'; 'whiteCat'; 'whiteCase'; 'noregretLearning'; 'optimization'}, 2, 'axis')    
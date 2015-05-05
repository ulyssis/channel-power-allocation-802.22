function convergenplot(sumUtilityWhitecat, sumUtilityWhitecase, sumUtilityNoregret)
hold off;    
figure(9);
    sumUtilityWhitecat_original = sumUtilityWhitecat;
    maxlength = max([size(sumUtilityWhitecat,2), size(sumUtilityWhitecase,2), size(sumUtilityNoregret,2)]);
    sumUtilityWhitecat(end+1:maxlength) = sumUtilityWhitecat(end);
    sumUtilityWhitecase(end+1:maxlength) = sumUtilityWhitecase(end);
    sumUtilityNoregret(end+1:maxlength) = sumUtilityNoregret(end);
    x = (1:maxlength);
    
    figure(30)
plot1 = plot (x, [sumUtilityWhitecat; sumUtilityWhitecase; sumUtilityNoregret]);
set(plot1(1),'Marker','.','DisplayName','whiteCat');
set(plot1(2),'Marker','+','DisplayName','whiteCase');
set(plot1(3),'Marker','d','DisplayName','noregret');
legend('whiteCat','whiteCase','NoregretLearning');            
legend('boxoff');
% title('convergence');
xlabel('Number of updates');
 ylabel('The sum of utility over all WBSs');
 
 figure(31)
plot2 =  plot ((1:size(sumUtilityWhitecat_original, 2)), sumUtilityWhitecat_original, '-b.');
%legend('whiteCat','whiteCase','NoregretLearning');


function [] = test_plot_group_errorbar(y, std_dev, num, numberOfchannelInUse)


c = 1:num;

%title('Title'); xlabel('x-label'); ylabel('y-label');
hold on
%%Bar(s)
%You can not color differently the same bar.
for i = 1:num
    bar(c(i)-0.15,y(i,1),0.2, 'b');
    bar(c(i),     y(i,2),0.2, 'r');
    bar(c(i)+0.15,y(i,3),0.2, 'g');
end
%%Errorbar
errH1 = errorbar(c-0.15,y(:,1),std_dev(:,1),'.','Color','k');
errH2 = errorbar(c,     y(:,2),std_dev(:,2),'.','Color','k');
errH3 = errorbar(c+0.15,y(:,3),std_dev(:,3),'.','Color','k');
errH1.LineWidth = 1.5;
errH2.LineWidth = 1.5;
errH3.LineWidth = 1.5;
%errH1.Color = [1 0.5 0];
%errH2.Color = [1 0.3 1];
%%Set x-ticks
xlim([c(1)-0.5, c(end) + 0.5]);


 %text(8, 1, strcat(num2str(numberOfchannelInUse), ' channels in use'));

function [] = plotMaxPowercdf(lp, cvx)
figure(11)
h1 = cdfplot(lp);
set(h1,'Color','r');
% set(h1,'Marker','+');
set(h1, 'LineStyle', '--');
hold on
h2 = cdfplot(cvx);
set(h2,'Color','b');
% set(h2,'Marker','*');
set(h2, 'LineStyle', '-');
axis([0 45 0 1]);

l = legend([h1, h2], 'Linear programming', 'convex programming');
set (l, 'Location', 'Southeast');

% legend(h1, 'Linear programming');
% hold all
% legend(h2, 'convex programming');

title('Cumulative distribution of maximal transmission power from different schemes');


    % Convert y-axis values to percentage values by multiplication
    a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
    % Create a vector of '%' signs
    pct = char(ones(size(a,1),1)*'%'); 
    % Append the '%' signs after the percentage values
    new_yticks = [char(a),pct];
    % 'Reflect the changes on the plot
    set(gca,'yticklabel',new_yticks);
    
    
    xlabel('Transmission power levels (W)');
%     ylabel('100% percentage');

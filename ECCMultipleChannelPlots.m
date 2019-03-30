


function ECCMultipleChannelPlots(plotLog, n, simSesult, runtimes, xstick)

% the latter snrRatio_random, snrRatio_whitecat, snrRatio_whitecase,
% snrRatio_noregret are sinr on each WBS in the last run, which is a
% intersection for sinr in one run

%% Average Tx
figure(plotLog+1);

numChannelScale = xstick(end) - xstick(1) + 1;
SchemesNum = 2;

% --- centralized Tx power
cenTx = simSesult(:, 1:n);
meanCenTx = zeros(1, numChannelScale);
stdCenTx = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(cenTx((i-1)*runtimes+1:i*runtimes, :), 2);
    meanCenTx(i) = mean(averageThR_ET);
    stdCenTx(i) = 1.96*std(averageThR_ET',1,2)/sqrt(runtimes);
end
% --- distributed Tx power
disTx = simSesult(:, n+1:2*n);
meanDisTx = zeros(1, numChannelScale);
stdDisTx = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(disTx((i-1)*runtimes+1:i*runtimes, :), 2);
    meanDisTx(i) = mean(averageThR_ET);
    stdDisTx(i) = 1.96*std(averageThR_ET',1,2)/sqrt(n);
end

% ---greedy Tx algorithm
greedyTx = simSesult(1:runtimes, 4*n+1 : 5*n);
averageTx_greedy_allRuns = mean(greedyTx, 2);
meanGreedyTx = mean(averageTx_greedy_allRuns);
stdGreedyTx = 1.96*std(averageTx_greedy_allRuns',1,2)/sqrt(n);

meanTx = [meanCenTx ; meanDisTx];
stdTx = [stdCenTx ; stdDisTx];

    h = gobjects(SchemesNum + 1, 1);
    for i = 1: SchemesNum
        x = xstick(1)+0.03*(i-1): 1 : xstick(end)+0.03*(i-1);
        y = meanTx(i, :);
        err = stdTx(i, :);
        h(i) = errorbar(x,y,err, '', '-o');
        hold on;
    end
    x = xstick(end) + 0.3;
    y = meanGreedyTx;
    err = stdGreedyTx;   
    h(SchemesNum + 1) = errorbar(x,y,err, '', 'gs', 'MarkerEdgeColor','green','MarkerFaceColor','green');
    
    LH = legend(h, {'Optimization', 'WhiteCat', 'Greedy Algo.'}, 'Location','northwest', 'FontSize', 16, 'Color', 'w', 'Box', 'on', 'EdgeColor', 'none');
    xticks(xstick);
    xlim([xstick(1)-0.2, xstick(end) + 0.5]);
    xlabel('Number of Channels for Transmission');
    ylabel('Avg. Tx Power per WBS (Watts)');
    set(gca,'FontSize',16);


    yl = ylim;
    lineH = line([xstick(end)+0.2 xstick(end)+0.2],yl,'Color',[0.9290, 0.6940, 0.1250], 'Linestyle', '--');
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    lineH.DisplayName = '';
    lineH.HandleVisibility = 'off';


    
    

    
%% capacity
figure(plotLog+2);

numChannelScale = size(simSesult, 1)/runtimes;
SchemesNum = 2;

% --- centralized capacity
cenCapacity = simSesult(:, 2*n+1:3*n);
meanCenCapacity = zeros(1, numChannelScale);
stdCenCapacity = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(cenCapacity((i-1)*runtimes+1:i*runtimes, :), 2);
    meanCenCapacity(i) = mean(averageThR_ET);
    stdCenCapacity(i) = 1.96*std(averageThR_ET',1,2)/sqrt(runtimes);
end

% --- distributed capacity
disCapacity = simSesult(:, 3*n+1:4*n);
meanDisCapacity = zeros(1, numChannelScale);
stdDisCapacity = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(disCapacity((i-1)*runtimes+1:i*runtimes, :), 2);
    meanDisCapacity(i) = mean(averageThR_ET);
    stdDisCapacity(i) = 1.96*std(averageThR_ET',1,2)/sqrt(n);
end

% ---greedy Tx algorithm
greedyCap = simSesult(1:runtimes, 5*n+1 : 6*n);
averageCap_greedy_allRuns = mean(greedyCap, 2);
meanGreedyCap = mean(averageCap_greedy_allRuns);
stdGreedyCap = 1.96*std(averageCap_greedy_allRuns',1,2)/sqrt(n);

meanCapacity = [meanCenCapacity ; meanDisCapacity];
stdCapacity = [stdCenCapacity ; stdDisCapacity];

    h = gobjects(SchemesNum + 1, 1);
    for i = 1: SchemesNum
        x = xstick(1)+0.03*(i-1): 1 : xstick(end)+0.03*(i-1);
        y = meanCapacity(i, :);
        err = stdCapacity(i, :);
        h(i) = errorbar(x,y,err, '', '-o');
        hold on;
    end
    x = xstick(end) + 0.3;
    y = meanGreedyCap;
    err = stdGreedyCap;   
    h(SchemesNum + 1) = errorbar(x,y,err, '', 'gs', 'MarkerEdgeColor','green','MarkerFaceColor','green');
    
    legend(h, {'Optimization', 'WhiteCat', 'Greedy Algo.'}, 'Location','northwest', 'FontSize', 16, 'Color', 'w', 'Box', 'on', 'EdgeColor', 'none');
    xticks(xstick);
    xlim([xstick(1)-0.2, xstick(end) + 0.5]);
    xlabel('Number of Channels for Transmission');
    ylabel('Throughput (bits/s)');
    set(gca,'FontSize',16);
    
    yl = ylim;
    lineH = line([xstick(end)+0.2 xstick(end)+0.2],yl,'Color',[0.9290, 0.6940, 0.1250], 'Linestyle', '--');
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    lineH.DisplayName = '';
    lineH.HandleVisibility = 'off';

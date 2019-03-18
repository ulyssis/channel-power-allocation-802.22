


function ECCMultipleChannelPlots(plotLog, n, simSesult, runtimes, xstick)

% the latter snrRatio_random, snrRatio_whitecat, snrRatio_whitecase,
% snrRatio_noregret are sinr on each WBS in the last run, which is a
% intersection for sinr in one run

%% Average Tx
figure(plotLog+1);

numChannelScale = xstick(end) - xstick(1) + 1;
SchemesNum = 2;

cenTx = simSesult(:, 1:n);
meanCenTx = zeros(1, numChannelScale);
stdCenTx = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(cenTx((i-1)*runtimes+1:i*runtimes, :), 2);
    meanCenTx(i) = mean(averageThR_ET);
    stdCenTx(i) = 1.96*std(averageThR_ET',1,2)/sqrt(runtimes);
end

disTx = simSesult(:, n+1:2*n);
meanDisTx = zeros(1, numChannelScale);
stdDisTx = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(disTx((i-1)*runtimes+1:i*runtimes, :), 2);
    meanDisTx(i) = mean(averageThR_ET);
    stdDisTx(i) = 1.96*std(averageThR_ET',1,2)/sqrt(n);
end

meanTx = [meanCenTx ; meanDisTx];
stdTx = [stdCenTx ; stdDisTx];

    h = gobjects(SchemesNum, 1);
    for i = 1: SchemesNum
        x = xstick(1)+0.03*(i-1): 1 : xstick(end)+0.03*(i-1);
        y = meanTx(i, :);
        err = stdTx(i, :);
        h(i) = errorbar(x,y,err);
        hold on;
    end
    legend(h, {'Optimization', 'WhiteCat'}, 'Location','northwest', 'FontSize', 16, 'Color', 'w', 'Box', 'on', 'EdgeColor', 'none');
    xticks(xstick);
    xlim([xstick(1)-0.2, xstick(end) + 0.2]);
    xlabel('Number of Channels for Transmission');
    ylabel('Avg. Tx Power per WBS (Watts)');
    set(gca,'FontSize',16);


    
    %% capacity
figure(plotLog+2);

numChannelScale = size(simSesult, 1)/runtimes;
SchemesNum = 2;

cenCapacity = simSesult(:, 2*n+1:3*n);
meanCenCapacity = zeros(1, numChannelScale);
stdCenCapacity = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(cenCapacity((i-1)*runtimes+1:i*runtimes, :), 2);
    meanCenCapacity(i) = mean(averageThR_ET);
    stdCenCapacity(i) = 1.96*std(averageThR_ET',1,2)/sqrt(runtimes);
end

disCapacity = simSesult(:, 3*n+1:4*n);
meanDisCapacity = zeros(1, numChannelScale);
stdDisCapacity = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(disCapacity((i-1)*runtimes+1:i*runtimes, :), 2);
    meanDisCapacity(i) = mean(averageThR_ET);
    stdDisCapacity(i) = 1.96*std(averageThR_ET',1,2)/sqrt(n);
end

meanCapacity = [meanCenCapacity ; meanDisCapacity];
stdCapacity = [stdCenCapacity ; stdDisCapacity];

    h = gobjects(SchemesNum, 1);
    for i = 1: SchemesNum
        x = xstick(1)+0.03*(i-1): 1 : xstick(end)+0.03*(i-1);
        y = meanCapacity(i, :);
        err = stdCapacity(i, :);
        h(i) = errorbar(x,y,err);
        hold on;
    end
    legend(h, {'Optimization', 'WhiteCat'}, 'Location','northwest', 'FontSize', 16, 'Color', 'w', 'Box', 'on', 'EdgeColor', 'none');
    xticks(xstick);
    xlim([xstick(1)-0.2, xstick(end) + 0.2]);
    xlabel('Number of Channels for Transmission');
    ylabel('Throughput (bits/s)');
    set(gca,'FontSize',16);
    


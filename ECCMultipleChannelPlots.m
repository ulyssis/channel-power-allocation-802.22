


function ECCMultipleChannelPlots(plotLog, n, simSesult, runtimes, xstick)

% the latter snrRatio_random, snrRatio_whitecat, snrRatio_whitecase,
% snrRatio_noregret are sinr on each WBS in the last run, which is a
% intersection for sinr in one run

%% Average Tx
figure(plotLog+1);

numChannelScale = xstick(end) - xstick(1) + 1;
SchemesNum = 3;

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
greedyTx = simSesult(:, 4*n+1 : 5*n);
meangreedyTx = zeros(1, numChannelScale);
stdgreedyTx = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(greedyTx((i-1)*runtimes+1:i*runtimes, :), 2);
    meangreedyTx(i) = mean(averageThR_ET);
    stdgreedyTx(i) = 1.96*std(averageThR_ET',1,2)/sqrt(n);
end
%---

meanTx = [meanCenTx ; meanDisTx; meangreedyTx];
stdTx = [stdCenTx ; stdDisTx; stdgreedyTx];
%---
    h = gobjects(SchemesNum, 1);
    for i = 1: SchemesNum
        x = xstick(1)+0.03*(i-1): 1 : xstick(end)+0.03*(i-1);
        y = meanTx(i, :);
        err = stdTx(i, :);
        h(i) = errorbar(x,y,err, '', '-o');
        hold on;
    end

    
    LH = legend(h, {'Optimization', 'WhiteCat', 'Greedy Algo.'}, 'Location','northwest', 'FontSize', 16, 'Color', 'w', 'Box', 'on', 'EdgeColor', 'none');
    xticks(xstick);
    xlim([xstick(1)-0.2, xstick(end) + 0.5]);
    xlabel('Number of Channels for Transmission');
    ylabel('Avg. Tx Power per WBS (Watts)');
    set(gca,'FontSize',16);


%     yl = ylim;
%     lineH = line([xstick(end)+0.2 xstick(end)+0.2],yl,'Color',[0.9290, 0.6940, 0.1250], 'Linestyle', '--');
%     ax = gca;
%     ax.XGrid = 'off';
%     ax.YGrid = 'on';
%     lineH.DisplayName = '';
%     lineH.HandleVisibility = 'off';


    
    

    
%% capacity
figure(plotLog+2);

numChannelScale = size(simSesult, 1)/runtimes;
SchemesNum = 3;

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
greedyCapacity = simSesult(:, 5*n+1 : 6*n);
meangreedyCapacity = zeros(1, numChannelScale);
stdgreedyCapacity = zeros(1, numChannelScale);
for i = 1:numChannelScale
    averageThR_ET = mean(greedyCapacity((i-1)*runtimes+1:i*runtimes, :), 2);
    meangreedyCapacity(i) = mean(averageThR_ET);
    stdgreedyCapacity(i) = 1.96*std(averageThR_ET',1,2)/sqrt(n);
end
%---

meanCapacity = [meanCenCapacity ; meanDisCapacity; meangreedyCapacity];
stdCapacity = [stdCenCapacity ; stdDisCapacity; stdgreedyCapacity];

    h = gobjects(SchemesNum, 1);
    for i = 1: SchemesNum
        x = xstick(1)+0.03*(i-1): 1 : xstick(end)+0.03*(i-1);
        y = meanCapacity(i, :);
        err = stdCapacity(i, :);
        h(i) = errorbar(x,y,err, '', '-o');
        hold on;
    end

    
    legend(h, {'Optimization', 'WhiteCat', 'Greedy Algo.'}, 'Location','northwest', 'FontSize', 16, 'Color', 'w', 'Box', 'on', 'EdgeColor', 'none');
    xticks(xstick);
    xlim([xstick(1)-0.2, xstick(end) + 0.5]);
    xlabel('Number of Channels for Transmission');
    ylabel('Avg. Shannon Capacity on End Terminals (bits/s)');
    set(gca,'FontSize',16);



%% Average Tx ratios of every WBS
figure(plotLog + 3);

numChannelScale = xstick(end) - xstick(1) + 1;
SchemesNum = 3;
averageThR_ET_allWBSs = zeros(numChannelScale, n); 
std_ET_allWBSs = zeros(numChannelScale, n); 
metaDataTxMean = zeros(n, (SchemesNum)*numChannelScale); % (scheme1, s2, s3)_w1, (scheme1, s2, s3)_w2 ...
metaDataTxStd = zeros(n, (SchemesNum)*numChannelScale);

% --- centralized Tx power
cenTx = simSesult(:, 1:n);

for i = 1:numChannelScale
    averageThR_ET_allWBSs(i, :) = mean(cenTx((i-1)*runtimes+1:i*runtimes, :), 1);
    std_ET_allWBSs(i, :) = (1.96*std(cenTx((i-1)*runtimes+1:i*runtimes, :),1,1)/sqrt(runtimes));
    metaDataTxMean(:, (i-1)*(SchemesNum)+1) = averageThR_ET_allWBSs(i, :)';
    metaDataTxStd(:, (i-1)*(SchemesNum)+1) = std_ET_allWBSs(i, :)';
end

% --- distributed Tx power
disTx = simSesult(:, n+1:2*n);

for i = 1:numChannelScale
    averageThR_ET_allWBSs(i, :) = mean(disTx((i-1)*runtimes+1:i*runtimes, :), 1);
    std_ET_allWBSs(i, :) = (1.96*std(disTx((i-1)*runtimes+1:i*runtimes, :),1,1)/sqrt(runtimes));
    metaDataTxMean(:, (i-1)*(SchemesNum)+2) = averageThR_ET_allWBSs(i, :)';
    metaDataTxStd(:, (i-1)*(SchemesNum)+2) = std_ET_allWBSs(i, :)';
end


% ---greedy Tx algorithm
greedyTx = simSesult(:, 4*n+1 : 5*n);

for i = 1:numChannelScale
    averageThR_ET_allWBSs(i, :) = mean(greedyTx((i-1)*runtimes+1:i*runtimes, :), 1);
    std_ET_allWBSs(i, :) = (1.96*std(greedyTx((i-1)*runtimes+1:i*runtimes, :),1,1)/sqrt(runtimes));
    metaDataTxMean(:, (i-1)*(SchemesNum)+3) = averageThR_ET_allWBSs(i, :)';
    metaDataTxStd(:, (i-1)*(SchemesNum)+3) = std_ET_allWBSs(i, :)';
end



h = linkaxis(metaDataTxMean, metaDataTxStd, n, SchemesNum); % numChannelScale = 4
set([h(1).XLabel],'string','Index of WBSs');
set([h(3).YLabel],'string','Avg. Tx Power per WBS (Watts)');

    
    
    %% capacities of every WBS
figure(plotLog + 4);



numChannelScale = xstick(end) - xstick(1) + 1;
SchemesNum = 3;
averageThR_ET_allWBSs = zeros(numChannelScale, n); 
std_ET_allWBSs = zeros(numChannelScale, n); 
metaDataTxMean = zeros(n, (SchemesNum)*numChannelScale); % (scheme1, s2, s3)_w1, (scheme1, s2, s3)_w2 ...
metaDataTxStd = zeros(n, (SchemesNum)*numChannelScale);

% --- centralized capacities
cenTx = simSesult(:, 2*n+1:3*n);

for i = 1:numChannelScale
    averageThR_ET_allWBSs(i, :) = mean(cenTx((i-1)*runtimes+1:i*runtimes, :), 1);
    std_ET_allWBSs(i, :) = (1.96*std(cenTx((i-1)*runtimes+1:i*runtimes, :),1,1)/sqrt(runtimes));
    metaDataTxMean(:, (i-1)*(SchemesNum)+1) = averageThR_ET_allWBSs(i, :)';
    metaDataTxStd(:, (i-1)*(SchemesNum)+1) = std_ET_allWBSs(i, :)';
end

% --- distributed capacities
disTx = simSesult(:, 3*n+1:4*n);

for i = 1:numChannelScale
    averageThR_ET_allWBSs(i, :) = mean(disTx((i-1)*runtimes+1:i*runtimes, :), 1);
    std_ET_allWBSs(i, :) = (1.96*std(disTx((i-1)*runtimes+1:i*runtimes, :),1,1)/sqrt(runtimes));
    metaDataTxMean(:, (i-1)*(SchemesNum)+2) = averageThR_ET_allWBSs(i, :)';
    metaDataTxStd(:, (i-1)*(SchemesNum)+2) = std_ET_allWBSs(i, :)';
end


% ---greedy capacities
greedyTx = simSesult(:, 5*n+1 : 6*n);

for i = 1:numChannelScale
    averageThR_ET_allWBSs(i, :) = mean(greedyTx((i-1)*runtimes+1:i*runtimes, :), 1);
    std_ET_allWBSs(i, :) = (1.96*std(greedyTx((i-1)*runtimes+1:i*runtimes, :),1,1)/sqrt(runtimes));
    metaDataTxMean(:, (i-1)*(SchemesNum)+3) = averageThR_ET_allWBSs(i, :)';
    metaDataTxStd(:, (i-1)*(SchemesNum)+3) = std_ET_allWBSs(i, :)';
end



h = linkaxis(metaDataTxMean, metaDataTxStd, n, SchemesNum); % numChannelScale = 4
set([h(1).XLabel],'string','Index of WBSs');
set([h(3).YLabel],'string', 'Avg. Shannon Capacity on End Terminals (bits/s)');
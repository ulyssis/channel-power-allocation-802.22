function [h] = linkaxis(metaDataTxMean, metaDataTxStd, n, SchemesNum)

figure
ax1 = subplot(4,1,1);
y = metaDataTxMean(:, 1:SchemesNum);
std_dev = metaDataTxStd(:, 1:SchemesNum);
num = n; %number of different subcategories
test_plot_group_errorbar(y, std_dev, num, 1);
title('one channel in use');


ax2 = subplot(4,1,2);
y = metaDataTxMean(:, SchemesNum+1:2*SchemesNum);
std_dev = metaDataTxStd(:, SchemesNum+1: 2*SchemesNum);
num = n; %number of different subcategories
test_plot_group_errorbar(y, std_dev, num, 2);
title('2 channels in use');

ax3 = subplot(4,1,3);
y = metaDataTxMean(:, 2*SchemesNum+1: 3*SchemesNum);
std_dev = metaDataTxStd(:, 2*SchemesNum+1: 3*SchemesNum);
num = n; %number of different subcategories
test_plot_group_errorbar(y, std_dev, num, 3);
title('3 channels in use');

ax4 = subplot(4,1,4);
y = metaDataTxMean(:, 3*SchemesNum+1: 4*SchemesNum);
std_dev = metaDataTxStd(:, 3*SchemesNum+1: 4*SchemesNum);
num = n; %number of different subcategories
test_plot_group_errorbar(y, std_dev, num, 4);
title('4 channels in use');

linkaxes([ax1,ax2,ax3,ax4],'x');
h = findobj(gcf,'type','axes');

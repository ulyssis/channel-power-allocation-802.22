data = [39.63,186.19;2.18,19.46];

b = bar(data);

ch = get(b,'children');

set(gca,'XTickLabel',{'LZCK[19]','Our scheme'})

legend('Application provider(AP)','Participant');

ylabel('Authentication time consumption(ms)');

 applyhatch(gcf,'|-/');

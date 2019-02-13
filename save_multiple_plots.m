FolderName = '/Users/max/Documents/git_li/channel-power-allocation-802.22';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
appendix = 0;
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, FigName, char("ECC_16_4_5e1-8_5000_1_20_10_" + appendix +".fig")));
  appendix = appendix + 1;
end
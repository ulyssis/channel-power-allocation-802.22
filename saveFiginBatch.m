% save multiple plots in the same time!
h=findobj('type','figure'); % find the handles of the opened figures
folder='/home/li/work/dev/DiCAPS/2011OtcGameOptReform/codes_cameraready/figures/';  % Desination folder
for k=1:numel(h)
  filename=sprintf('%d.pdf',h(k));
  file=fullfile(folder,filename);
  saveas(h(k),file);
end
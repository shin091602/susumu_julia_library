% reshape figure and save fig for paper and ppt
function [] = save_fig(fig,filename,view_angle)
% fig        : figure object
% filename   : filename excluding .eps
% view_angle : = [az,el]

  % change the view point
  ax = findobj(fig,'Type','axes');
  for j = 1:max(size(ax))
    view(ax(j),view_angle(1),view_angle(2));
  end

  % save
  ppt_name = strcat(filename,'.png');
  tunefig('document',fig);
  print(fig,'-dpng',ppt_name,'-r300');

  ppt_name = strcat(filename,'_ppt.png');
  tunefig('ppt',fig);
  print(fig,'-dpng',ppt_name,'-r300');

  characters = strlength(filename);
  if characters > 180
    % too long for windows
    throw(MException("main:save_fig","Error: filename is too long"));
  end

  doc_name = strcat(filename,'.pdf');
  tunefig('document',fig);
  exportgraphics(fig,doc_name,'ContentType','vector');
end
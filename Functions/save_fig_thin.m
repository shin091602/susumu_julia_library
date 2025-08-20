% reshape figure and save fig for paper and ppt
function [] = save_fig_thin(fig,filename,view_angle)
% fig        : figure object
% filename   : filename excluding .eps
% view_angle : = [azimuth, elevation]

  % change the view point
  ax = findobj(fig,'Type','axes');
  for j = 1:max(size(ax))
    view(ax(j),view_angle(1),view_angle(2));
  end

  % save
  ppt_name = strcat(filename,'.png');
  tunefig('thin_document',fig);
  print(fig,'-dpng',ppt_name,'-r300');

  ppt_name = strcat(filename,'_ppt.png');
  tunefig('thin_ppt',fig);
  print(fig,'-dpng',ppt_name,'-r300');

  characters = strlength(filename);
  if characters > 180
    % too long for windows
    throw(MException("main:save_fig","Error: filename is too long"));
  end

  doc_name = strcat(filename,'.pdf');
  tunefig('thin_document',fig);
  exportgraphics(fig,doc_name,'ContentType','vector');
end
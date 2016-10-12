function h=mifit_Help_About(fig)
% Help/About: display the About dialogue. The handle ID is in adddata(gcf, 'handle_About')
  if nargin ==0, fig=''; end
  icon = fullfile(ifitpath,'Docs','images','ILL-web-jpeg.jpg');
  
  % Display About dialog
  t = [ sprintf('Welcome to miFit, a GUI to iFit.\n ') version(iData,2) sprintf('.\n Visit <http://ifit.mccode.org>') ];
  if isempty(dir(icon))
    h = msgbox(t,'miFit: About','help');
  else
    h = msgbox(t,'miFit: About','custom', imread(icon));
  end
  g=findobj(h, 'type','uicontrol');
  config = getappdata(mifit_fig, 'Preferences');
  set(g,'fontsize',config.FontSize);
  if ~isempty(fig)
    setappdata(fig, 'handle_About', h);
  end

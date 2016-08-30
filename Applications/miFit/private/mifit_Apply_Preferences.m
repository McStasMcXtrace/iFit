function mifit_Apply_Preferences
  fig = mifit_fig;
  config = getappdata(fig, 'Preferences');
  % change all Font Size
  h=[ findobj(fig, 'Type','uicontrol') ; ...
      findobj(fig, 'Type','axes') ; findobj(fig, 'Type','text') ; ...
      findobj(fig, 'Type','uipanel') findobj(fig, 'Type','uitable') ];
  set(h, 'FontSize', config.FontSize);
  % for uimenu, could we use tip given by Y Altman 
  % <https://fr.mathworks.com/matlabcentral/newsreader/view_thread/148095>
  % h = findobj(fig, 'Type','uimenu');

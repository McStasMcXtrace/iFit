function ud=iData_plot_contextmenu(a, h, xlab, ylab, zlab,  T, S, d, cmd, mp, mname)
% iData_plot_contextmenu: add a contextmenu to the iData plot
% used in iData/plot

% ==============================================================================
% assemble UserData stuff

ud = help(a); % gather properties and other information 

if isempty(h), return; end
tproperties = ud.tproperties;
mproperties = ud.mproperties;

% get back any previous axis UserData
ud0 = get(gca, 'UserData');
if ~isempty(ud0) && isfield(ud0, 'handles')
  handles = ud0.handles;
else handles = [];
end

if iscell(handles), handles = [ handles{:} ]; end
handles = [ handles(:)' h(:)' ];

% treat HG2 hggroup which do not properly transfer UIContextMenu settings

% ==============================================================================
% contextual menu for the single object being displayed
% internal functions must be avoided as it uses LOTS of memory

try
    index1 = findobj(get(h,'Children'),'Tag',[ 'iData_plot_' a.Tag '_contextmenu_object' ]);
catch
    index1 = [];
end
index = [ findobj(h,'Tag',[ 'iData_plot_' a.Tag '_contextmenu_object' ]) ...
          findobj(gcf,'Tag',[ 'iData_plot_' a.Tag '_contextmenu_object' ]) ];
try
    index1 = findobj(get(h,'Children'),'Tag',[ 'iData_plot_' a.Tag '_contextmenu_object' ]);
    index = [ index index1 ];
end

% check if a context menu is already attached
if all(isempty(index))
  uicm = uicontextmenu('Tag',[ 'iData_plot_' a.Tag '_contextmenu_object' ]); 
  % menu About
  uimenu(uicm, 'Label', [ 'About ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ' ...' ], ...
    'Callback', [ 'msgbox(getfield(get(get(gco,''UIContextMenu''),''UserData''),''properties''),' ...
                  '''About: Figure ' num2str(double(gcf)) ' ' T ' <' S '>'',' ...
                  '''custom'',getfield(getframe(gcf),''cdata''), get(gcf,''Colormap''));' ] );

  % menu Toggle error bars (1D object)
  if ndims(a) == 1
    uimenu(uicm, 'Label','Toggle Error Bars', 'Callback', [...
     'tmp_h=gco; if strcmpi(get(get(tmp_h(1),''Parent''),''Type''), ''Hggroup''), ' ...
     'tmp_h=get(tmp_h,''Parent''); end; ' ...
     'tmp_h=findobj(tmp_h,''Type'',''errorbar''); ' ...
     'if isempty(tmp_h), tmp_h=get(gco,''children''); tmp_h=tmp_h(2:end); end; '...
     'if ~isempty(tmp_h), if strcmp(get(tmp_h(1),''visible''),''off''), tmp_v=''on''; else tmp_v=''off''; end;' ...
     'set(tmp_h,''visible'',tmp_v); end; clear tmp_h tmp_v' ]);
  end
  uimenu(uicm, 'Separator','on', 'Label', [ 'Title: "' T '" ' d ]);
  if exist(a.Source,'file') && ~isdir(a.Source)
    uimenu(uicm, 'Label', [ 'Source: <' S '>' ], 'Callback',[ 'edit(''' a.Source ''')' ]);
  else
    uimenu(uicm, 'Label', [ 'Source: <' S '>' ]);
  end
  if ~isempty(d)
    uimenu(uicm, 'Label', [ 'Label: ' d ]);
  end
  % menu List of Commands (history)
  uimenu(uicm, 'Label', [ 'Cmd: ' cmd ' ...' ], 'Callback', [ ...
   'tmp_ud = get(get(gco,''UIContextMenu''),''UserData'');' ...
   'listdlg(''ListString'', getfield(tmp_ud, ''commands''),' ...
    '''ListSize'',[400 300],''Name'', tmp_ud.title ,''PromptString'', tmp_ud.name );clear tmp_ud' ])
  uimenu(uicm, 'Label', [ 'User: ' a.User ]);

  uimenu(uicm, 'Separator','on', 'Label', '[Rank]         [Value] [Description]');
  for index=1:numel(tproperties)
    uimenu(uicm, 'Label', tproperties{index});
  end
  % model parameters
  if numel(mproperties) > 1
    for index=1:min(numel(mproperties),10)
      if index==1
      uimenu(uicm, 'Separator','on', 'Label', mproperties{index});
      else
      uimenu(uicm, 'Label', mproperties{index});
      end
    end
    if ~isempty(mproperties) && index < numel(mproperties)
      uimenu(uicm, 'Label', '...');
    end
  end
    
  % menu About iFit
  uimenu(uicm, 'Separator','on','Label', 'About iFit/iData', 'Callback', ...
    [ 'msgbox(''' version(iData,2) sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);

  ud.contextmenu_object = uicm;
  ud.handles = handles;
  
  set(uicm,'UserData', ud);
  set(gca, 'UserData', ud);
  set(h,   'UIContextMenu', uicm);  % flag the presence of the context menu
  
  % set the menu also inside group objects
  try
    if strcmp(get(h,'Type'),'hggroup')
      h = get(h,'Children');
    end
  catch
    h = [];
  end
  for index=1:numel(h)
    set(h(index),   'UIContextMenu', uicm);
  end
end

% ==============================================================================
% contextual menu for the axis frame
if ~isempty(get(gca, 'UserData'))
  ud = get(gca, 'UserData');
end
if isempty(get(gca,   'UIContextMenu'))
  uicm = uicontextmenu('Tag','iData_plot_contextmenu_gca');
  % menu Duplicate (axis frame/window)
  uimenu(uicm, 'Label', 'Duplicate View...', 'Callback', ...
     [ 'tmp_cb.g=gca; tmp_cb.ud=get(gca,''UserData'');' ...
       'tmp_cb.f=figure; tmp_cb.c=copyobj(tmp_cb.g,gcf); ' ...
       'set(tmp_cb.c,''position'',[ 0.1 0.1 0.85 0.8]);' ...
       'set(gcf,''Name'',''Copy of ' strrep(char(a),'''','"') '''); ' ...
       'set(gca,''XTickLabelMode'',''auto'',''XTickMode'',''auto'');' ...
       'set(gca,''YTickLabelMode'',''auto'',''YTickMode'',''auto'');' ...
       'set(gca,''ZTickLabelMode'',''auto'',''ZTickMode'',''auto'');' ...
       'title(tmp_cb.ud.title);', ...
       'xlabel(tmp_cb.ud.xlabel);ylabel(tmp_cb.ud.ylabel); clear tmp_cb;']);
       
  if ndims(a) == 1 && ~isfield(ud,'contextual_1d')
    ud.contextual_1d = 1;
  end
  % menu Toggle all error bars (axis)
  % we handle 'old' (Matlab < 2014b) errorbars (hggroup) by hidding the 2nd line
  % we handle tne 'new' Matlab >= 2014b by toggling between the errorbar and the single line
  if isfield(ud,'contextual_1d') && ud.contextual_1d==1
    uimenu(uicm, 'Label','Toggle All Error Bars', 'Callback', [ ... 
      'tmp_h = findobj(gca,''type'',''errorbar''); ' ...
      'if isempty(tmp_h), ' ...
      'tmp_hg=findobj(gca,''Type'',''hggroup''); ' ...
      'for tmp_i=1:numel(tmp_hg);' ...
      'tmp_v=get(tmp_hg(tmp_i),''Children''); ' ...
      'if numel(tmp_v) > 1, tmp_h = [ tmp_h tmp_v(2) ]; end; ' ...
      'end; end; tmp_v=[]; '...
      'if strcmp(get(tmp_h(1),''visible''),''off''), tmp_v=''on''; else tmp_v=''off''; end; ' ...
      'set(tmp_h,''visible'',tmp_v); ' ...
      'clear tmp_h tmp_v tmp_hg;' ...
     ]);
  end
  uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
  uimenu(uicm, 'Label','Show legend', 'Callback', ...
    [ 'tmp_cb.g=gca; tmp_cb.ud=get(gca,''UserData'');' ...
      'if isfield(tmp_cb.ud,''handles''), legend(tmp_cb.ud.handles,''location'',''best''); end;' ...
      'clear tmp_cb;' ] ...
  );
  if ndims(a) >= 2 && ~isfield(ud,'contextual_2d')
    ud.contextual_2d = 1;
  end
  % menus entries for lin/log and appearence
  if isfield(ud,'contextual_2d') && ud.contextual_2d==1
    uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback', [ ...
      '[tmp_a,tmp_e]=view; if (tmp_a==0 & tmp_e==90) view(3); else view(2); end;' ...
      'clear tmp_a tmp_e; lighting none;alpha(1);shading flat;rotate3d off;axis tight;' ]);
    uimenu(uicm, 'Label','Smooth View','Callback', 'shading interp;');
    uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
    uimenu(uicm, 'Label','Add Transparency','Callback', 'alphamap(''decrease''); for tmp_h=get(gca, ''children'')''; try; alpha(tmp_h,0.7*get(tmp_h, ''facealpha'')); end; end; clear tmp_h');
    uimenu(uicm, 'Label',[ 'Linear/Log signal ' strtok(title(a)) ],...
      'Callback', 'if strcmp(get(gca,''zscale''),''linear'')  set(gca,''zscale'',''log''); else set(gca,''zscale'',''linear''); end');
    uimenu(uicm, 'Label',[ 'Linear/Log X axis ' strtok(xlabel(a)) ], ...
      'Callback', 'if strcmp(get(gca,''xscale''),''linear'')  set(gca,''xscale'',''log''); else set(gca,''xscale'',''linear''); end');
    uimenu(uicm, 'Label',[ 'Linear/Log Y axis ' strtok(ylabel(a)) ], ...
      'Callback', 'if strcmp(get(gca,''yscale''),''linear'')  set(gca,''yscale'',''log''); else set(gca,''yscale'',''linear''); end');
    uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
  else
    uimenu(uicm, 'Label','Reset View', 'Callback','view(2);lighting none;alpha(1);shading flat;axis tight;rotate3d off;');
    uimenu(uicm, 'Label',[ 'Linear/Log signal ' strtok(title(a)) ],'Callback', 'if strcmp(get(gca,''yscale''),''linear'')  set(gca,''yscale'',''log''); else set(gca,''yscale'',''linear''); end');
    uimenu(uicm, 'Label',[ 'Linear/Log axis ' strtok(xlabel(a)) ],'Callback', 'if strcmp(get(gca,''xscale''),''linear'')  set(gca,''xscale'',''log''); else set(gca,''xscale'',''linear''); end');
  end

  uimenu(uicm, 'Separator','on','Label', 'About iFit/iData', ...
    'Callback',[ 'msgbox(''' version(iData,2) sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);

  ud.contextmenu_gca = uicm;
  set(gca, 'UserData',      ud);
  set(gca, 'UIContextMenu', uicm);
end

% ==============================================================================
% contextual menu for the figure
% add rotate/pan/zoom tools to the figure in case java machine is not started
if ~usejava('jvm')
  if isempty(get(gcf, 'UIContextMenu'))
    uicmf = uicontextmenu('Tag','iData_plot_contextmenu_fig');
    uimenu(uicmf, 'Label','Zoom on/off', 'Callback','zoom');
    uimenu(uicmf, 'Label','Pan on/off',  'Callback','pan');
    if ndims(a) >= 2
      uimenu(uicmf, 'Label', 'Rotate on/off', 'Callback','rotate3d');
    end
    uimenu(uicmf, 'Label','Legend on/off', 'Callback','legend(gca, ''toggle'',''Location'',''Best'');');
    uimenu(uicmf, 'Label','Print...', 'Callback','printpreview');

    set(gcf, 'UIContextMenu', uicmf);
    set(gcf, 'KeyPressFcn', @(src,evnt) eval('if lower(evnt.Character)==''r'', lighting none;alpha(1);shading flat;axis tight;rotate3d off; zoom off; pan off; end') );
    ud.contextmenu_fig = uicmf;
  end
end



function popupmessage(filename)
% popupmessage: basic text editor written in 100% Matlab

% the Figure
handles(1)=figure('units','pixels',...
    'position',[250 250 700 700],...
    'menubar','none');

handles(2)=uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[0.1 0.01 0.1 0.05],...
    'string','Load ...', 'ForegroundColor','blue',...    
    'callback',@event_load);

% the Text/Edit
handles(3)=uicontrol('style','edit',...
    'units','normalized',...
    'position',[0.01 0.07 0.98 0.92],...
    'string','',...
    'HorizontalAlignment','left',...
    'max',300);
    
handles(4)=uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[0.25 0.01 0.1 0.05],...
    'string','Save as...', 'ForegroundColor','blue',...    
    'callback',@event_save);
    
handles(5)=uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[0.85 0.01 0.1 0.05],...
    'string','Close', 'ForegroundColor','red',...  
    'Tooltip','Close this window without saving content.', ...  
    'callback',@event_close);
    
handles(6)=uicontrol('style','text',...
    'units','normalized',...
    'position',[0.4 0.01 0.4 0.05],...
    'string','');
  
if nargin == 0
  filename = '';
end
if ~isempty(filename)
  action_load(filename);
end
  
% ----------------------------------------------------------------------------

  function event_load(obj,event)
      action_load('');
  end
  
  function action_load(filename)
    if nargin == 0,       filename == ''; end
    if isempty(filename), filename = uigetfile; end
    if ~ischar(filename) || all(filename == 0),     return; end
    
    if ~isempty(dir(filename))
      content=fileread(filename);
      titl = filename;
    else
      content = filename;
      titl = content(:)';
    end
    set(handles(3),'string',content);
    
    if length(titl) > 80, titl = [ titl(1:79) ' ...' ]; end
    if ~isempty(titl)
      % update information text in lower part
      set(handles(1), 'name', titl);
      set(handles(6), 'ToolTip', [ 'File:' titl sprintf('\nSize:') num2str(length(content)) ], 'String',titl);
    end
  end
  
  function event_save(obj,event)
      action_save('');
  end
  
  function action_save(filename)
    if nargin == 0,       filename == ''; end
    if isempty(filename), filename = uiputfile; end
    if ~ischar(filename) || all(filename == 0),     return; end

    content = get(handles(3),'string');
    fid = fopen(filename);
    if fid == -1
      error([ mfilename ': Could not open file ' filename ]);
    end
    fprintf(fid, '%s', content);
    fclose(fid);
    % update information text in lower part
    set(handles(1), 'name', filename);
    set(handles(6), 'ToolTip', [ 'File:' filename sprintf('\nSize:') num2str(length(content)) ], 'String',titl);
  end
  
  function event_close(obj,event)
      delete(handles(1));
  end
  
end



function fallback_web(url)
% fallback_web: basic text editor written in 100% Matlab

  list = {};

  if nargin == 0
    url = '';
  end
  if isempty(url), url = [ ifitpath filesep 'Docs' filesep 'index.html' ]; end

  % use Java browser
  if usejava('jvm') 
      if ~isempty(dir(url))
        url = [ 'file://' url ]; % local file
      end
      root       = fileparts(url);
      if strncmp(root, 'file://',7)
        root = root(8:end);
      end
      handles(1) = figure('menubar','none', 'Name', url, 'KeyPressFcn',@KeyPressFcn);
      % button bar [ Back URL History Home ]
      % Back button
      handles(2) = uicontrol('style','pushbutton','Units','Normalized', ...
        'ToolTip', 'Go back one page', ...
        'String','<-','Position',[0.01 0.91 0.04 0.05], 'Callback',@setBack);
      % URL text
      handles(3) = uicontrol('style','edit','Units','Normalized', ...
        'ToolTip', 'Type a URL here', 'BackgroundColor','white', ...
        'String', url, 'Position',[0.06 0.91 0.63 0.05], 'Callback',@setURL);
      % History
      handles(4) = uicontrol('style','pushbutton','Units','Normalized', ...
        'ToolTip', 'List of previous URL', ...
        'String', 'History', 'Position',[0.70 0.91 0.14 0.05], 'Callback',@setFromHistory);
      % HOME button
      handles(5) = uicontrol('style','pushbutton','Units','Normalized', ...
        'ToolTip', sprintf('Go back to iFit main page\n%s', [ ifitpath filesep 'Docs' filesep 'index.html' ] ), ...
        'String','HOME','Position',[0.85 0.91 0.14 0.05], 'Callback',@setHome);
        
      je         = javax.swing.JEditorPane('text/html',url);
      jp         = javax.swing.JScrollPane(je);
      [hcomponent, hcontainer] = javacomponent(jp, [], handles(1));
      set(hcontainer, 'units', 'normalized', 'position', [0,0,1,.9]);
      
      je.setEditable(false);
      setPage(url);
      set(je, 'HyperlinkUpdateCallback',@action_follow_link)
  else
    disp('Can not display Web page (Java not available). Open it manually :')
    disp(url)
  end

% ------------------------------------------------------------------------------
%                              Navitation actions
% ------------------------------------------------------------------------------
  function action_follow_link(obj,event)
    % triggered when passing over a link in the JEditorPane
    l    = get(obj);                      % this is a structure
    event= l.HyperlinkUpdateCallbackData; % this is the event caught by the Listener
    if strcmp(event.eventType,'ACTIVATED')
      link = event.description;
      setPage(link);
    end
  end
    
  function setPage(link)
    % used to set the URL in the JEditorPane
    if any(link == '#') % is there an chor here ?
      [link, anchor] = strtok(link, '#');
      if isempty(anchor), link = [ '#' link ]; end
    else
      anchor = '';
    end
    if ~isempty(dir([ root filesep link ]))
      link = [ root filesep link ];
    elseif link(1) == '#'
      link = [ url link ];
    end
    if ~isempty(dir(link))
      link = [ 'file://' link ]; % local file
    end
    je.setPage([ link anchor ]);
    url        = strtok(link,'#');
    root       = fileparts(url);
    if strncmp(root, 'file://',7)
      root = root(8:end);
    end
    list{end+1} = url;
  end
  
  function KeyPressFcn(src,evnt)
    % used to handle key pressed in the figure
    % 'h'         -> setHome
    % 'backspace' -> back in history
    if lower(evnt.Key) == 'h'
      setHome
    elseif lower(evnt.Key) == 'backspace'
      
    end
  end
  
  function setHome(src,evnt)
    setPage([ ifitpath filesep 'Docs' filesep 'index.html' ]);
  end
  
  function setURL(src, evnt)
    link = get(handles(3), 'String');
    if ~isempty(link)
      setPage(link);
    else
      set(handles(3), 'String', url);
    end
  end
  
  function setFromHistory(src, evnt)
    link = listdlg('ListString',list, 'InitialValue', length(list), ...
      'ListSize', [ 400 100 ], 'SelectionMode', 'single', ...
      'Name','Select URL to go back to');
    if isempty(link), return; end
    setPage(list{link})
    list = list(1:link);
  end
  
  function setBack(src, evnt)
    if length(list) > 1
      l = length(list);
      setPage(list{l-1});
      list = list(1:(l-1));
    end
  end
end

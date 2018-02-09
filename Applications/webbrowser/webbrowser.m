function ret=webbrowser(url, method)
% webbrowser: basic web browser running from Matlab
%
%   Opens a simplistic web browser, built from Matlab/Java.
%     Can be used as replacement for the 'web' command in deployed applications.
%     Requires a running JVM.
%   This browser has very limited rendering capabilities. It does not support 
%     JavaScript, and other 'modern' HTML extensions (flash, HTML5,...).
%   The browser has a display pane, a Home button, a Back button, an editable URL 
%   field, and keeps track of the navigation history. Does not support proxy settings.
%
%   webbrowser        Opens an empty browser.
%   webbrowser(file)  Display the specified file (full path or in the Matlab search path)
%   webbrowser(url)   Opens a browser displaying the specified URL.
%   webbrowser(url,'system') tries first to open the web page with the system browser.
%
%   Copyright: Licensed under the BSD, version 1.0.1
%              E. Farhi, ILL, France <farhi@ill.fr> Aug 2017, http://ifit.mccode.org

  list = {};

  if nargin == 0
    url = '';
  end
  if nargin < 2, method='matlab'; end
  
  if exist('ifitpath')
    Home = fullfile(ifitpath,'Docs','index.html');
  else
    Home = 'http://ifit.mccode.org';
  end
  
  if isempty(url), url = Home; end
  
  if isempty(dir(url)) 
    if ~isempty(which([ url '.html' ])) && feature('ShowFigureWindows') && usejava('jvm')
      url = which([ url '.html' ]); % prefer HTML over TXT when Java/Display
    elseif ~isempty(which([ url '.txt' ]))
      url = which([ url '.txt' ]); % prefer TXT when no Java/no Display
    elseif ~isempty(which([ url '.html' ]))
      url = which([ url '.html' ]); % use HTML if no TXT
    else
      url1 = which(url);
      if ~isempty(url1), url=url1; end
    end
  end
  
  handles=[];
  
  % return code from Browser: 0=OK, 1=error.
  ret=1;  % error code
  if strcmp(method, 'system')
    ret = open_system_browser(url);
  end
  
  % when fails, we open our own browser
  if ret && usejava('jvm') && usejava('swing') && feature('ShowFigureWindows')
      % use Java browser
      if ~isempty(dir(url)) % local file/dir
          if isdir(url), P=url; F=''; E='';
          else [P,F,E] = fileparts(url); end
          if isempty(P), P=pwd; end
          url = fullfile(P, [F E]);
          url = [ 'file://' url ]; % local file
      end
      root       = fileparts(url);
      if strncmp(root, 'file://',7)
        root = fileparts(url(8:end));
      end
      handles(1) = figure('menubar','none', 'Name', url,'NextPlot','new');
      % button bar [ Back URL History Home ]
      % Back button
      handles(2) = uicontrol('style','pushbutton','Units','Normalized', ...
        'ToolTip', 'Go back one page', ...
        'String','<-','Position',[0.01 0.91 0.04 0.05], 'Callback',@setBack);
      % URL text
      handles(3) = uicontrol('style','edit','Units','Normalized', ...
        'ToolTip', 'Type a URL here (http://, file://, ...) or "matlab:<command>"', 'BackgroundColor','white', ...
        'String', url, 'Position',[0.06 0.91 0.63 0.05], 'Callback',@setURL);
      % History
      handles(4) = uicontrol('style','pushbutton','Units','Normalized', ...
        'ToolTip', 'List of previous URL', ...
        'String', 'History', 'Position',[0.70 0.91 0.14 0.05], 'Callback',@setFromHistory);
      % HOME button
      handles(5) = uicontrol('style','pushbutton','Units','Normalized', ...
        'ToolTip', sprintf('Go back to Home main page\n%s',Home ), ...
        'String','HOME','Position',[0.85 0.91 0.14 0.05], 'Callback',@setHome);
      % File Menu
      handles(6) = uimenu(handles(1), 'Label', 'File');
      uimenu(handles(6), 'Label', 'Open URL', 'Callback', @openURL, 'Accelerator','O');
      uimenu(handles(6), 'Label', 'Save as...', 'Callback', @savePage, 'Accelerator','S');
      uimenu(handles(6), 'Label', 'Export page', 'Callback', @printPage, 'Accelerator','P');
      uimenu(handles(6), 'Separator','on', 'Label', 'Quit', ...
        'Callback', @exitBrowser, 'Accelerator','Q');
      % Navigation Menu
      handles(6) = uimenu(handles(1), 'Label', 'Navigation');
      uimenu(handles(6), 'Label', 'Home', 'Callback', @setHome, 'Accelerator','H');
      uimenu(handles(6), 'Label', 'Back', 'Callback', @setBack, 'Accelerator','B');
      uimenu(handles(6), 'Label', 'Set from history...', 'Callback', @setFromHistory, 'Accelerator','L');
      uimenu(handles(6), 'Separator','on', 'Label', 'About this Browser', 'Callback', @aboutBrowser);
      
      je         = javax.swing.JEditorPane('text/html',url);
      jp         = javax.swing.JScrollPane(je);
      [hcomponent, hcontainer] = javacomponent(jp, [], handles(1));
      set(hcontainer, 'units', 'normalized', 'position', [0,0,1,.9]);
      je.setDragEnabled(true)
      je.setEditable(false);
      ret=setPage(url);
      set(je, 'HyperlinkUpdateCallback',@action_follow_link);
      set(handles(1), 'UserData', je);
  end
  if ret % could not open browser
    disp(url)
    if ~isempty(dir(url)), type(url); end
  end
  
  function ret=open_system_browser(url)
  % opens URL with system browser. Returns non zero in case of error.
    if strncmp(url, 'file://', length('file://'))
      url = url(8:end);
    end
    ret = 1;
    if ismac,  precmd = 'DYLD_LIBRARY_PATH= ;';
    elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
    else precmd=''; end
    if ispc
      ret=system([ precmd 'start "' url '"']);
    elseif ismac
      ret=system([ precmd 'open "' url '"']);
    else
      [ret, message]=system([ precmd 'xdg-open "' url '"']);
    end
  end

% ------------------------------------------------------------------------------
%                              Navigation actions
% ------------------------------------------------------------------------------
  function action_follow_link(obj,event)
    % triggered when passing over a link in the JEditorPane
    l    = get(obj);                      % this is a structure
    event= l.HyperlinkUpdateCallbackData; % this is the event caught by the Listener
    if isfield(event, 'eventType') && strcmp(event.eventType,'ACTIVATED')
      link = event.description;
      setPage(link);
    end
  end
  
  function exitBrowser(src,evnt) % Menu File:Quit
    if ~isempty(handles)
      delete(handles(1));
    end
  end
  
  function printPage(src,evnt)   % Menu File:Print
    filemenufcn(gcbf,'FileSaveAs');
  end
  
  function savePage(src,evnt)    % Menu File:Save as
    filename = uiputfile('*.html','Save page as HTML');
    if ~ischar(filename) || all(filename == 0),     return; end
    content  = je.getText;
    fid = fopen(filename, 'w+');
    if fid == -1
      error([ mfilename ': Could not open file ' filename ]);
    end
    fprintf(fid, '%s', char(content));
    fclose(fid);
  end
    
  function ret=setPage(link)         % set the JEditorPane content (URL)
    % ret: return code from Browser: 0=OK, 1=error.
    ret=0;
    % used to set the URL in the JEditorPane
    if any(link == '#') % is there an anchor here ?
      [link, anchor] = strtok(link, '#');
      if isempty(anchor), link = [ '#' link ]; end
    else
      anchor = '';
    end
    % test if this is a 'matlab:' command
    if strncmp(link, 'matlab:', length('matlab:'))
      url = link; root = link;
      % check if we have %20 in the URL and replace them by spaces
      if ~isempty(strfind(link, '%20'))
        link = strrep(link, '%20', ' ');
      end
      if strncmp(link, 'matlab:helpwin', length('matlab:helpwin')) ...
        || strncmp(link, 'matlab:doc', length('matlab:doc'))
        link = [ 'matlab:web ' link(15:end) ]; % uses embedded help system
      end
      % evaluate the matlab command
      try
        link = evalc(link(8:end));
      end
      link = strrep(link, sprintf('\n'),'<br>');
      % highlight the first line
      index = strfind(link, '<br>');
      if length(index) > 1
        index=index(1);
        link = [ '<b>' link(1:(index-1)) ' </b> ' link(index:end) ];
      end
      % put the result of the matlab command in the browser
      try
        je.setText( link );
      end
    else % this is a supported URL
      url = get(handles(3), 'String');  % the current URL
      if ~isempty(dir([ root filesep link ]))
        link = [ root filesep link ];
      elseif link(1) == '#'
        link = [ url link ];
      end
      if ~isempty(dir(link))
        link = [ 'file://' link ]; % local file
      end
      % check if we have spaces in the URL and replace them by %20
      if any(link == ' ')
        link = strrep(link, ' ', '%20');
      end
      try
        je.setPage([ link anchor ]);
      catch ME % usually java.net.MalformedURLException from JEditorPane
        ret=1;
      end
      % in case of error, try with http:// in front if missing
      if ret
        try
          je.setPage([ 'http://' link anchor ]);
          link = [ 'http://' link anchor ];
          ret=0;
        catch ME
          ret=1;
        end
      end
        
      % in case of error, try with the initial root.
      if ret
        try
          je.setPage([ root '/' link anchor ]);
          link = [ root '/' link anchor ];
          ret=0;
        catch ME
          ret=1;
        end
      end
      % in case of error, try with the initial url as relative root.
      if ret
        try
          je.setPage([ url '/' link anchor ]);
          link = [ url '/' link anchor ];
          ret=0;
        catch ME
          ret=1;
        end
      end
      if ret
        disp([ mfilename ': Can not open URL ' link ])
        disp(getReport(ME))
        return
      end
      url        = strtok(link,'#');
      root       = fileparts(url);
      % root should only contain the dns
      split = textscan(url, '%s','Delimiter','/'); 
      if ~isempty(split)
        split=split{1};
        split=split(cellfun(@isempty, split));
        if ~isempty(split) && any(strcmp(split{1}, {'http:','https','file:','ftp:'}))
          split=split(2:end);
        end
        if ~isempty(split)
          root = split{1};
        end
      end
      if strncmp(root, 'file://',7)
        root = root(8:end);
      end
    end
    list{end+1} = url;
    set(handles(3), 'String', url);
  end
  
  function setHome(src,evnt)     % go back to Home
    setPage(Home);
  end
  
  function openURL(src,evnt)     % open a URL from a Dialog
    link = inputdlg('Enter the URL to go to','Open URL',1,{''});
    if ~isempty(link)
      setPage(link{1});
    end
  end
  
  function setURL(src, evnt)     % set URL by changing the URL edit box content
    if ~isempty(handles)
      link = get(handles(3), 'String');
      if ~isempty(link)
        setPage(link);
      else
        set(handles(3), 'String', url); % restore the previous URL
      end
    end
  end
  
  function setFromHistory(src, evnt)  % list the history and jump to item
    link = listdlg('ListString',list, 'InitialValue', length(list), ...
      'ListSize', [ 400 100 ], 'SelectionMode', 'single', ...
      'Name','Select URL to go back to');
    if isempty(link), return; end
    setPage(list{link})
    list = list(1:link);
  end
  
  function setBack(src, evnt)    % go back one URL
    if length(list) > 1
      l = length(list);
      setPage(list{l-1});
      list = list(1:(l-1));
    end
  end
  
  function aboutBrowser(src, evnt) % about dialog
    helpdlg(sprintf('This is a simplistic web browser, built from Matlab/Java. Does not support proxy settings. E. Farhi, ILL, France <farhi@ill.fr> Aug 2017. Copyright: Licensed under BSD. VISIT http://ifit.mccode.org for more.'), 'About this simplistic Browser');
  end
  
  % protect figure from over-plotting
  if ~isempty(handles)
    set(handles(1),'HandleVisibility','callback');
  end
  
end

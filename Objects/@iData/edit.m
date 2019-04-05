function h = edit(a, option)
% edit(s) : edit iData object Signal/Monitor
%
%   @iData/edit function to view the Signal/Monitor of a data set. 
%     The data appears as a spreadsheet (requires Java to be enabled).
%
%   An additional 'option' argument allows to specify what to do with the Table
%   d=edit(iData, handle)   The given Table handle is exported into an iData object.
%   edit(iData, 'editable') The Table data can be modified after creation.
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=edit(a);
%
% Version: $Date$ $Version$ $Author$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

% the Selection is stored in the figure UserData
  h = [];
  if nargin < 2
    option = '';
  end
  
  if ishandle(option)
    h = iData_edit_iData(option);
    return
  end

  if numel(a) > 1
    for index=1:numel(a)
      h(index) = edit(a(index), option);
    end
    return
  end
  
  if ndims(a) > 2
    warning([ mfilename ': iData object ' a.Tag ' dimension is ' num2str(ndims(a)) ' but should be 1 or 2 to be edited. Skipping.' ]);
    h = 0;
    return
  end

  if exist('uitable') && isjava('jvm')
    if prod(size(a)) > 1e5
      iData_private_warning(mfilename, [ 'Object ' a.Tag ' is too large (numel=' num2str(prod(size(a))) ...
    '.\n\tYou should rebin with e.g. a=a(1:2:end, 1:2:end, ...).' ]);
    end
    
    if strfind(option, 'editable')
      editable = true;
    else
      editable = false;
    end
    
    Signal    = getaxis(a, 'Signal'); % Signal/Monitor
    
    % opens a Table with the Signal content there-in
    NL = sprintf('\n');
    f = figure('Name', [ 'Edit ' char(a) ], 'MenuBar','none'); % raw Signal, no monitor weightening
    p = get(f, 'Position');
    h = uitable('Data', Signal);
    
    
    tooltip = sprintf([ 'Here is the Data set, without Axes definition.\n' ...
      'Access the contextual menu with the mouse right-click:\n' ...
      '* You can select a subset and plot it.\n' ...
      '* You can export its content into a new iData object, or in the clipboard.']);
    if editable
      tooltip = [ tooltip sprintf('\n* You can modify the Table content.\n* You can paste the clipboard (matrix or file name).') ];
    end
    
    % UserData.Monitor = real(get(a,'Monitor'));
    
    UserData.Selection = [];
    UserData.handle    = h;
    UserData.figure    = f;
    UserData.Name      = char(a);
    UserData.Axes      = getaxis(a, []);
    UserData.modified  = false;
    
    set(h, 'Position', [0,0,p(3),p(4)]); 

    set(h, 'Tag',[ mfilename '_' a.Tag ]);
    set(h, 'Units','normalized');
    set(h, 'Position', [0,0,1,1]);
    set(h, 'CellSelectionCallback', @iData_edit_CellSelectionCallback);
    set(h, 'ToolTipString', tooltip);
    if editable, set(h, 'ColumnEditable', true); end

    set(f, 'UserData', UserData);  % contains the selection indices
    set(f, 'NextPlot', 'new');
    
    % set the Axes as ColumnName and RowName
    if ndims(a) == 1
      % ColumnName == getaxis(a,1)
      x = unique(getaxis(a,1));
      if numel(x) == numel(Signal)
        x = cellfun(@(c){num2str(c)}, num2cell(x));
        x{1} = strtrim([ xlabel(a) ' ' x{1} ]);
        if size(a,2) == 1
          set(h, 'RowName', x);
        else
          set(h, 'ColumnName', x);
        end
      end
    elseif ndims(a) == 2
      % ColumnsName == getaxis(a,2)
      y = unique(getaxis(a,2));
      if numel(y) == size(a,2)
        y=cellfun(@(c){num2str(c)}, num2cell(y));
        y{1} = strtrim([ xlabel(a) ' ' y{1} ]);
        set(h, 'ColumnName', y);
      end
      % RowName     == getaxis(a,1)
      x = unique(getaxis(a,1));
      if numel(x) == size(a,1)
        x = cellfun(@(c){num2str(c)}, num2cell(x));
        x{1} = strtrim([ ylabel(a) ' ' x{1} ]);
        set(h, 'RowName', x);
      end
    end

%      'CellEditCallback', @iData_edit_CellEditCallback, , ...
%      'RearrangeableColumn','on', 'ColumnEditable', true, ...
%      'ToolTipString', [ char(a) NL 'You may edit the Signal from the object.' ...
%         NL 'The Contextual menu enables to plot a selection of the oject' ], ...    
    
    % attach contextual menu
    uicm = uicontextmenu; 
    S=a.Source; T=a.Title;
    if length(T) > 23, T=[ T(1:20) '...' ]; end
    if length(S) > 23, S=[ '...' S(end-20:end) ]; end
    uimenu(uicm, 'Label', [ 'Title: "' T '"' ]);
    uimenu(uicm, 'Label', [ 'Source: <' S '>' ]);
    uimenu(uicm, 'Label', [ 'Tag: '  a.Tag  ]);
    uimenu(uicm, 'Label', [ 'User: ' a.User ]);
    if ~isempty(a.Label)
      uimenu(uicm, 'Label', [ 'Label: ' a.Label ]);
    end
    if length(getaxis(a))
      uimenu(uicm, 'Separator','on', 'Label', '[Rank]         [Value] [Description]');
      for index=0:length(getaxis(a))
        [v, l] = getaxis(a, num2str(index));
        x      = getaxis(a, index);
        m      = get(a, 'Monitor');
        if index==0 & not(all(m==1 | m==0))
          uimenu(uicm, 'Label', sprintf('%6i %15s  %s [%g:%g] (per monitor)', index, v, l, min(x(:)), max(x(:))));
        else
          uimenu(uicm, 'Label', sprintf('%6i %15s  %s [%g:%g]', index, v, l, min(x(:)), max(x(:))));
        end
      end
    end
    uimenu(uicm, 'Separator','on','Label','Copy selection to clipboard', 'Callback', @iData_edit_copy);
    if editable
    uimenu(uicm, 'Label','Paste here from clipboard', 'Callback', @iData_edit_paste);
    uimenu(uicm, 'Label','Resize data set', 'Callback', @iData_edit_resize);
    end
    uimenu(uicm, 'Label','Export to main workspace...', 'Callback', @iData_edit_export);
    uimenu(uicm, 'Label','Export to CSV file...', 'Callback', @iData_edit_export_csv);
    uimenu(uicm, 'Label','Plot selection...', 'Callback', @iData_edit_display);
    uimenu(uicm, 'Separator','on','Label', 'About iData', 'Callback',[ 'msgbox(''' version(iData) ''')' ]);
    % attach contexual menu to the table
    try
      set(h,   'UIContextMenu', uicm); 
    end
    
  else
    % open the data file
    try
      h = a.Source;
      edit(a.Source);
    end
  end


% ------------------------------------------------------------------------------
% Events
% ------------------------------------------------------------------------------

function [Selection, rows, columns, UserData] = iData_edit_getselection
  % get the current selection in the uitable
  fig = gcbf;
  UserData  = get(fig, 'UserData');
  Signal    = get(UserData.handle, 'Data');
  Selection = UserData.Selection; % as { rows, columns }
  
  % get the Selection as Signal( rows, columns )
  if isempty(Selection) || isempty(Selection{1}) || isempty(Selection{2})
    Selection = Signal;
    rows    = [ 1 size(Signal,1) ];
    columns = [ 1 size(Signal,2) ];
  else
    rows    = Selection{1};
    columns = Selection{2};
    Selection = Signal(Selection{:});
  end
  if isempty(rows)   || any(rows    <= 0), rows=1; end
  if isempty(columns)|| any(columns <= 0), columns=1; end

% function called when the Table selection is plotted
function iData_edit_display(hObject, eventdata, varargin)
  % hObject is the handle of the context menu item
  
  [Selection, rows, columns] = iData_edit_getselection;
  UserData  = get(fig, 'UserData');
  % plot it (no axes here ?)
  if isscalar(Selection) || ndims(Selection) > 2
    return;
  end
  nfig = figure;
  if isvector(Selection)
    plot(Selection);
    if numel(rows) == 1
      xlabel(sprintf('(%i,[%i:%i])', rows, min(columns),max(columns)));
    else
      xlabel(sprintf('([%i:%i],%i)', min(rows),max(rows), columns));
    end
  elseif ndims(Selection) == 2
    surf(Selection);
    ylabel(sprintf('[%i:%i]', min(rows),max(rows)));
    xlabel(sprintf('[%i:%i]', min(columns),max(columns)));
  end
  t = sprintf('Edit ([%i:%i,%i:%i]) %s', ...
        min(rows),max(rows),min(columns),max(columns), UserData.Name);
  set(nfig, 'Name', t);
  title(t)
  

function iData_edit_CellSelectionCallback(source, eventdata)
  % copy the selection into the UserData
  
  fig = gcbf; % current UItable figure
  
  rows    = unique(eventdata.Indices(:,1)); % rows and columns selected
  columns = unique(eventdata.Indices(:,2));
  
  UserData = get(fig, 'UserData');
  UserData.Selection = { rows,columns };
  set(fig, 'UserData',UserData);

% Clipboard management
function iData_edit_copy(hObject, eventdata, varargin)
  % copy the selection to the clipboard
  
  [Selection, rows, columns, UserData] = iData_edit_getselection;
  
  if isempty(Selection), return; end
  
  disp([ mfilename ': copy Data Selection ' mat2str(size(Selection)) ' to the clipboard.' ]);
  disp(UserData.Name);
  c = num2str(Selection);
  c = cellstr(c);
  c = sprintf('%s\n', c{:});
  
  clipboard('copy', c);
  
function iData_edit_paste(hObject, eventdata, varargin)
  % paste from the clipboard at the current location
  
  [Selection, rows, columns, UserData] = iData_edit_getselection;
  
  Selection = clipboard('paste');
  if isempty(Selection), return; end
  
  if ~isempty(dir(Selection)) % a file name
    Selection = iData(Selection);
    Selection = getaxis(Selection, 'Signal');
  else
  
    % split into lines
    Selection = regexp(Selection,'\n|\r','split');
    Selection = Selection(~cellfun(@isempty, Selection));
    Selection = str2num(str2mat(Selection));  % now a double array
  end
  
  if isempty(Selection), return; end
  
  fig = gcbf;
  Data      = get(UserData.handle, 'Data');
  
  % test if the selection will replace some data without extending borders
  sz_Data = size(Data);
  sz_Sel  = size(Selection);
  sz_Ins  = [ min(rows) min(columns) ]; % insertion location
  rows    = sz_Ins(1):(sz_Sel(1)+sz_Ins(1)-1);
  columns = sz_Ins(2):(sz_Sel(2)+sz_Ins(2)-1);
  if any(sz_Data < sz_Sel+sz_Ins-1)
    % the Selection will extend the Data set
    
    newData = zeros(sz_Sel+sz_Ins-1);
    % put the Data in
    newData(1:sz_Data(1),1:sz_Data(2)) = Data;
    Data = newData;
    clear newData
    % this also invalidates the stored Axes from initial object
    UserData.Axes = [];
    set(fig, 'UserData', UserData)
  end
  Data( rows, columns ) = Selection;
  set(UserData.handle, 'Data', Data);
  UserData.modified = true;
  set(fig, 'UserData', UserData);
  
function iData_edit_export(varargin)
  % exports the data set to the 'ans' variable in main workspace
  fig = gcbf;
  UserData  = get(fig, 'UserData');
  d=iData_edit_iData(UserData.handle);
  if isempty(d), return; end
  disp([ mfilename ': export Data ' mat2str(size(d)) ' into iData object ''ans''.' ]);
  disp(UserData.Name);

  assignin('base','ans',d);
  % make 'ans' be displayed
  ans = d

function iData_edit_export_csv(varargin)
  % exports the data set to a CSV file
  fig = gcbf;
  UserData  = get(fig, 'UserData');
  d=get(UserData.handle, 'Data');
  if isempty(d), return; end
  filename = uiputfile(...
    {'*.csv','Comma Separated Value file (Excel, *.csv)';'*.*','All files'}, ...
    [ mfilename ': Save as CSV' ]);
  if isempty(filename), return; end
  csvwrite(filename, double(d));

function d=iData_edit_iData(h)
  % create an iData object from Table
  d = [];
  if ~ishandle(h) || ~strcmpi(get(h, 'Type'),'uitable')
    return
  end
  
  % get the Data
  fig       = get(h,'Parent');
  UserData  = get(fig, 'UserData');
  if isfield(UserData,'modified') && ~UserData.modified
    return
  end
  Signal    = get(h, 'Data');
  
  % get the Axes
  Axes      = UserData.Axes;
  % check axes if there is at least one given
  if ~isempty(Axes)
    d = iData(Axes{:},Signal);
  else
    d = iData(Signal);
  end
  d.Title = [ 'edit(' UserData.Name ')' ];
  
function iData_edit_resize(varargin)
  % resize the Data set, padding with zeros
  
  fig = gcbf;
  UserData  = get(fig, 'UserData');
  h = UserData.handle;
  
  % get current Data
  Data    = get(h, 'Data');
  sz_Data = size(Data);
  % display Dialog
  prompt = {'New number of rows:','New number of columns:'};
  name   = [ mfilename ': Expand/reduce the Table' ];
  options.Resize='on';
  options.WindowStyle='normal';
  defaultanswer = {num2str(sz_Data(1)), num2str(sz_Data(2)) };
  answer=inputdlg(prompt,name,1,defaultanswer);
  
  if numel(answer) < 2, return; end
  sz_nData=[ str2num(answer{1}) str2num(answer{2}) ];
  if numel(sz_nData) ~= 2, return; end
  
  newData = zeros(sz_nData);
  rows    = min([ sz_nData(1) sz_Data(1) ]);
  columns = min([ sz_nData(2) sz_Data(2) ]);
  newData(rows, columns) = Data(rows, columns);
  set(h, 'Data', newData);
  UserData.modified = true;
  set(fig, 'UserData', UserData);


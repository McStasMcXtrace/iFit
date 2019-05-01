function uitable(s,varargin)
% UITABLE A dialogue which allows to modify objects in a table
%   UITABLE(s) edits a single object or array. The property values can be changed
%   and are updated when closing the window. The table has a contextual menu.
%   New properties can be added. The table can be reverted to its initial 
%   content, and printed. The dialog window is 'modal' by default.
%   A modal dialog box prevents a user from interacting with other windows before
%   responding to the modal dialog box.
%
%   UITABLE(s, OPTIONS) specifies additional options for the table layout and
%   behaviour. The argument OPTIONS is a structure which can contain fields
%     Name            The name of the dialogue (string)
%     FontSize        The FontSize used to display the table (default: system FontSize).
%     ListString      The labels to be used for each structure field (cell).
%                     The cell must have items 'fieldname description...' 
%                     (default: use the property names).
%     TooltipString   A string to display (as help)
%     CreateMode      Can be 'modal' (wait, default) or 'non-modal' (return immediately).
%                     When non-modal, objects are updated when closing the window.
%     Tag             A tag for the dialogue.
%     ColumnFormat    A cell which specifies how to handle the property values.
%                     'char' all values can be anything.
%                     'auto' protects scalar numeric values from being changed to something else (default).
%
%   UITABLE(s, 'option1', value1, ...) is the same as above, but specified options
%   as separate name/value pairs.
%
% In modal (default) creation mode, the updated object(s) is/are returned upon
%   closing the window. A modal dialog box prevents a user from interacting
%   with other windows before responding to the modal dialog box.
% In non-modal creation mode, the window is displayed, and remains visible
%   while the execution is resumed. The call to uitable returns a
%   configuration which should be used as follows to retrieve the modified
%   structure later:
%     ad = uitable(a, struct('CreateMode','non-modal'));
%     % continue execution, and actually edit and close the window
%     % ...
%     % get the modified structure (when closing window)
%     Data = getappdata(0,ad.tmp_storage);
%     a    = cell2struct(reshape(Data(:,2:end), ad.size), Data(:,1), 1);
%     rmappdata(0, ad.tmp_storage);
%
% input:
%   structure: the initial struct to edit
%   options:   a set of options, namely:

%
% Example:
%   a.Test=1; a.Second='blah'; uitable(a)
%   options.ListString={'Test This is the test field','Second 2nd'};
%   uitable([a a], options);
%
% Version: $Date$ $Version$ $Author$

  fields  = fieldnames(s);    % members of the structure

  % get options or default values
  if nargin == 0, s=[]; end
  if isempty(s) || ~isstruct(s), return; end

  % get/build options
  if nargin ==2 && isstruct(varargin{1})
    options = varargin{1};
  elseif mod(nargin, 2) == 1 % name/value pairs
    options = [];
    for index=1:2:(nargin-1)
      if isvarname(varargin{index}) && index<nargin-1
        options.(varargin{index}) = varargin{index+1};
      end
    end
  end

  if ~isfield(options, 'Name')
    if ~isempty(inputname(1))
      options.Name = [ 'Edit ' class(s) ' ' inputname(1) ];
    else
      options.Name = [ 'Edit ' class(s) ' members' ];
    end
  end
  if ~isfield(options, 'FontSize')
    options.FontSize = get(0,'DefaultUicontrolFontSize');
  end
  if ~isfield(options, 'ListString')
    options.ListString = {};
  end
  if ~isfield(options, 'TooltipString')
    options.TooltipString = '';
  end
  if ~isfield(options, 'ColumnFormat')
    options.ColumnFormat = 'auto';
  end
  if ~isempty(options.TooltipString)
    options.TooltipString = [ options.TooltipString sprintf('\n') ];
  end
  options.TooltipString = [ options.TooltipString ...
    'The configuration will be updated when you close this window.' sprintf('\n') ...
     'You can edit the Field names and subsequent values.' sprintf('\n') ...
     'To Cancel edition, use the CANCEL context menu item (right click).' ];
  [~,tmp_storage] = fileparts(tempname);
  if ~isfield(options,'Tag')
    options.Tag = [ mfilename '_' tmp_storage ];
  end
  if ~isfield(options,'CreateMode')
    options.CreateMode = 'modal';
  end

  % create a uitable from the structure fields and values
  f = figure('Name',options.Name, 'MenuBar','none', ...
      'NextPlot','new', ...
      'Tag', options.Tag, 'Units','pixels', ...
      'CloseRequestFcn', @uitable_close);

  % override default mechanism for Table update when closing
  if isfield(options, 'CloseRequestFcn') && ~strcmp(options.CreateMode, 'modal')
    set(f, 'CloseRequestFcn', options.CloseRequestFcn);
  end

  % sort ListString so that it matches the fields
  ListString = cell(size(fields));
  for index=1:numel(fields)
    [tokens, rems] = strtok(options.ListString);
    rems = strtrim(rems);
    index_f = find(strcmp(fields{index}, tokens),1);
    if ~isempty(index_f) && ~isempty(rems{index_f})
      ListString{index} = rems{index_f};
    else
      ListString{index} = fields{index};
    end
  end
  options.ListString = ListString;
  structure = rmfield(struct(s), 'Private'); % convert to struct

  % create the Table content to display, handle array of structures
  % check if the structure values are numeric, logical, or char
  Data0       = cell(numel(fields), numel(structure)+2);  % 1st column=ListString, % 2nd column=fieldnames
  fields_type = cell(size(Data0));
  Data0(:,1)  = ListString(:);
  Data0(:,2)  = fields(:);
  
  for index_s = 1:numel(structure)
    for index=1:numel(fields)
      Data0{index, index_s+2} = getfield(structure(index_s), fields{index});
      item = Data0{index,index_s+2};
      % uitable only support char or scalar numeric/logical
      flag = ischar(item) || (numel(item) ==1 && (isnumeric(item) || islogical(item)));
      if ~flag && exist('class2str')
        fields_type{index,index_s+2} = class(item);
        Data0{index,index_s+2} = class2str('', Data0{index,index_s+2}, 'eval');
      elseif strcmp(options.ColumnFormat,'char')
        fields_type{index,index_s+2} = class(item);
        if ~ischar(item) && ~isobject(item) && isnumeric(item)
          Data0{index,index_s+2} = num2str(item);
        end
      end
    end
  end

  % determine the window size to show
  % height is given by the number of fields
  height = (numel(options.ListString)+3)*options.FontSize*2;
  % width is given by the length of the longest RowName + nb of elements in array (columns)
  sz = cellfun(@numel,options.ListString);
  width  = (max(median(sz),mean(sz))+numel(structure)*5)*options.FontSize;
  % compare to current window size
  p = get(f, 'Position');
  p(3) = width;
  if p(4) > height, p(4) = height; end
  set(f, 'Position',p);
  % now we check that the width is not too small. 10 cm minimum.
  set(f, 'Units', 'centimeters');
  p = get(f, 'Position');
  if p(3) < 10, p(3) = 10; end
  set(f, 'Position',p);
  set(f, 'Units', 'pixels');

  % set ColumnName
  ColumnName = 'Description';
  ColumnName = [ ColumnName 'Field' num2cell(1:numel(structure)) ];
  ColumnEditable = [ false true ];
  ColumnWidth = {max(cellfun(@numel,options.ListString)+4)*options.FontSize/2, ...
    max(cellfun(@numel,fields)+4)*options.FontSize/2};
  ColumnFormat = { 'char' 'char' }; % Description and field names
  for index_s = 1:numel(structure)
    ColumnEditable(end+1) = true;
    ColumnWidth{end+1}    = 'auto';
    if strcmp(options.ColumnFormat,'char')
      ColumnFormat{end+1} = 'char';
    else
      ColumnFormat{end+1} = [];
    end
  end

  % create the table
  t = uitable('Parent',f, ...
    'Data',         Data0, ...
    'ColumnEditable',ColumnEditable, ...
    'Tag',          [ options.Tag '_Table' ], ...
    'FontSize',     options.FontSize, ...
    'Units','normalized', 'Position', [0.02 0.02 0.98 0.98 ], ...
    'ColumnWidth',  ColumnWidth, ...
    'TooltipString',options.TooltipString, ...
    'ColumnName',   ColumnName, ...
    'ColumnFormat', ColumnFormat);

  % add context menu to add a field
  uicm = uicontextmenu;
  uimenu(uicm, 'Label','Revert (undo all changes)', 'Callback', @uitable_revert);
  uimenu(uicm, 'Label','Append property/field',     'Callback', @uitable_append);
  uimenu(uicm, 'Label','OK',     'Callback', 'close(gcbf)', 'separator','on');
  uimenu(uicm, 'Label','CANCEL', 'Callback', 'delete(gcbf)');
  uimenu(uicm, 'Label','Print',  'Callback', 'printdlg(gcbf)');
  % attach contexual menu to the table
  try
    set(t,   'UIContextMenu', uicm);
  end

  % assemble the dialogue information structure, store UserData
  ad.figure       = f;
  ad.table        = t;
  ad.options      = options;
  ad.s0           = s;      % initial object(s)
  ad.Data0        = Data0;  % initial table content
  ad.fields       = fields; % this is also Data0(:,2)
  ad.fields_type  = fields_type;
  ad.size         = size(struct2cell(structure));   % initial Data size
  set(f, 'UserData', ad);

  % when the figure is deleted, we should get the uitable Data back
  % wait for figure close
  if strcmp(options.CreateMode, 'modal')
    uiwait(f);
  end

  % ==============================================================================
  function uitable_revert(source, evnt)
  % UITABLE_REVERT restore initial data
    ad = get(gcbf, 'UserData');
    set(ad.table, 'Data', ad.Data0);
  end

  function uitable_append(source, evnt)
  % UITABLE_APPEND add a new line to the table
    ad   = get(gcbf, 'UserData');
    Data = get(ad.table,'Data');
    fields_type = ad.fields_type;
    % append a line
    Data(end+1,:) = cell(1,size(Data,2));
    Data{end,1}   = 'New field';
    Data{end,2}   = 'New';
    fields_type(end+1,:) = cell(1,size(Data,2));
    % update storage
    ad.fields_type = fields_type;
    set(ad.table, 'Data', Data);
    set(gcbf, 'UserData', ad);
  end % uitable_append
  
  function uitable_close(source, evnt)
  % UITABLE_CLOSE update initial array of objects on exit
    ad   = get(gcbf, 'UserData');
    if ~isempty(ad)
      Data = get(ad.table,'Data');
      fields = Data(:,2); % field names
      for index_s = 1:numel(ad.s0)      % loop on objects
        for index=1:numel(fields)       % loop on fields
          val = Data{index, index_s+2};
          if ~isempty(fields_type{index, index_s+2})
            try
              val = eval(val);
            end
          end
          if numel(ad.s0) == 1, s = ad.s0;
          else                  s = ad.s0(index_s); end
          try
            set(s, fields{index}, val);
          end
        end
      end
    end
    delete(gcbf);
  end % uitable_close

end

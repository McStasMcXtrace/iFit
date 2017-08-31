function structure = structdlg(structure,options)
  % structdlg: a dialogue which allows to modify structure(s) in a table
  %
  % In modal (default) creation mode, the updated structure is returned upon
  %   closing the window. A modal dialog box prevents a user from interacting 
  %   with other windows before responding to the modal dialog box.
  % In non-modal creation mode, the window is displayed, and remains visible
  %   while the execution is resumed. The call to structdlg returns a
  %   configuration which should be used as follows to retrieve the modified 
  %   structure later:
  %     ad = structdlg(a, struct('CreateMode','non-modal'));
  %     % continue execution, and actually edit and close the window
  %     % ...
  %     % get the modified structure (when closing window)
  %     Data = getappdata(0,ad.tmp_storage); fields=Data(:,1); Data=Data(:,2:end);
  %     new_a      = cell2struct(reshape(Data, ad.size), fields, 1);
  %     rmappdata(0, ad.tmp_storage);
  %
  % input:
  %   structure: the initial struct to edit
  %   options:   a set of options, namely:
  %     options.Name: the name of the dialogue (string)
  %     options.FontSize: the FontSize used to display the table
  %       default: use the system FontSize.
  %     options.ListString: the labels to be used for each structure field (cell).
  %       the cell must have items 'fieldname description...'
  %       default: use the structure member names.
  %     options.TooltipString: a string to display (as help)
  %     options.CreateMode: can be 'modal' (default) or 'non-modal'.
  %     options.Tag: a tag for the dialogue.
  %     options.CloseRequestFcn: a function handle or expresion to execute when 
  %       closing the dialogue in 'non-modal' mode. 
  %
  % output:
  %   structure: the modified structure in 'modal' mode (default),
  %     or the dialogue information structure in 'non-modal' mode.
  %
  % Example: 
  %   a.Test=1; a.Second='blah'; structdlg(a)
  %   options.ListString={'Test This is the test field','Second 2nd'};
  %   structdlg([a a], options);
  %
  % Version: $Date$
  % (c) E.Farhi, ILL. License: BSD. <ifit.mccode.org>
  
  fields  = fieldnames(structure);    % members of the structure
  
  % get options or default values
  if nargin == 0, structure=[]; end
  if isempty(structure) || ~isstruct(structure), return; end
  if nargin == 1, options=[]; end
  
  
  if ~isfield(options, 'Name')
    if ~isempty(inputname(1))
      options.Name = [ 'Edit ' inputname(1) ];
    else
      options.Name = [ 'Edit structure members' ];
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
  if ~isempty(options.TooltipString)
    options.TooltipString = [ options.TooltipString sprintf('\n') ];
  end
  options.TooltipString = [ options.TooltipString ...
    'The configuration will be updated when you close this window.' ];
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
      'CloseRequestFcn', ...
        [ 'setappdata(0,''' tmp_storage ''', get(findobj(gcbf,''Tag'',''' options.Tag '_Table' '''),''Data'')); delete(gcbf)' ]);
  
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
  
  % create the Table content to display, handle array of structures
  % check if the structure values are numeric, logical, or char
  Data0       = cell(numel(fields), numel(structure)+2);  % 1st column=ListString, % 2nd column=fieldnames
  fields_type = cell(size(Data0));
  Data0(:,1)  = ListString(:);
  Data0(:,2)  = fields(:);
  for index_s = 1:numel(structure)
    Data0(:,index_s+2) = struct2cell(structure(index_s));
    for index=1:numel(fields)
      item = Data0{index,index_s+2};
      % uitable only support char or scalar numeric/logical
      flag = ischar(item) || (numel(item) ==1 && (isnumeric(item) || islogical(item)));
      if ~flag && exist('class2str')
        fields_type{index,index_s+2} = class(item);
        Data0{index,index_s+2} = class2str('', Data0{index,index_s+2}, 'eval');
      end
    end
  end

  % determine the window size to show
  TextWidth = 12;
  TextHeight= 30;
  % height is given by the number of fields
  height = (numel(options.ListString)+3)*TextHeight;
  % width is given by the length of the longest RowName
  width = max(cellfun(@numel,options.ListString))*TextWidth + numel(structure)*options.FontSize*5;
  % compare to current window size
  p = get(f, 'Position');
  p(3) = width;
  if p(4) > height, p(4) = height; end
  set(f, 'Position',p);
  
  % set ColumnName
  ColumnName = 'Description';
  ColumnName = [ ColumnName 'Field' num2cell(1:numel(structure)) ];
  ColumnFormat='char';
  ColumnFormat=[ ColumnFormat 'char' cell(1,numel(structure)) ];
  ColumnEditable = [ false true ];
  ColumnWidth = {max(cellfun(@numel,options.ListString))*options.FontSize/2, ...
    max(cellfun(@numel,fields))*options.FontSize/2};
  for index_s = 1:numel(structure)
    ColumnEditable(end+1) = true;
    ColumnWidth{end+1}    = 'auto';
  end

  % create the table
  t = uitable('Parent',f, ...
    'Data', Data0, ...
    'ColumnEditable',ColumnEditable, ...
    'Tag', [ options.Tag '_Table' ], ...
    'FontSize',options.FontSize, ...
    'Units','normalized', 'Position', [0.05 0.2 .9 .7 ], ...
    'ColumnWidth',ColumnWidth,'TooltipString',options.TooltipString, ...
    'ColumnName',ColumnName);
    
  % add context menu to add a field
  uicm = uicontextmenu; 
  uimenu(uicm, 'Label','Revert (undo all changes)', 'Callback', @structdlg_revert);
  uimenu(uicm, 'Label','Append property/field',     'Callback', @structdlg_append);
  % attach contexual menu to the table
  try
    set(t,   'UIContextMenu', uicm); 
  end
  
  % assemble the dialogue information structure
  ad.figure       = f;
  ad.table        = t;
  ad.options      = options;
  ad.tmp_storage  = tmp_storage;
  ad.structure0   = structure;
  ad.fields       = fields;
  ad.size         = size(struct2cell(structure));   % initial Data
  set(f, 'UserData', ad);
  setappdata(0,tmp_storage,Data0);   % we shall collect the Table content here
  
  % when the figure is deleted, we should get the uitable Data back
  % wait for figure close
  if strcmp(options.CreateMode, 'modal')
    uiwait(f);
    Data = getappdata(0,tmp_storage);
    
    if isempty(Data), return; end
    % get new field names
    fields = genvarname(strtok(Data(:,2)));
    Data   = Data(:,3:end);
    fields_type = fields_type(:,3:end);
    % restore initial structure field type when set
    if prod(ad.size) == numel(Data)
      Data = reshape(Data, ad.size);
    end
    for index=1:numel(Data)
      if ~isempty(fields_type{index})
        try
          Data{index} = eval(Data{index});
        catch
          disp([ mfilename ': failed evaluation of ' fields{index} ' new value "' Data{index} '". Skipping.' ]);
        end
      end
    end

    % assemble new options...
    structure = cell2struct(Data, fields, 1);
    rmappdata(0, tmp_storage);
  else
    structure = ad;
  end
  
  % ==============================================================================
  function structdlg_revert(source, evnt)
  % restore initial data
    ad = get(gcbf, 'UserData');
    Data = getappdata(0,ad.tmp_storage);
    set(ad.table, 'Data',Data);
  end
    
  function structdlg_append(source, evnt)
  % add a new line to the table
    ad   = get(gcbf, 'UserData');
    Data = get(ad.table,'Data');
    % append a line
    Data(end+1,:) = cell(1,size(Data,2));
    fields_type(end+1,:) = cell(1,size(Data,2));
    Data{end,1}   = 'New field'; 
    Data{end,2}   = 'New'; 
    set(ad.table, 'Data',Data);
  end
  
end

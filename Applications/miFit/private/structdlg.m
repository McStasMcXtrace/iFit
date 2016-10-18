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
  %     % get the modified structure
  %     new_a      = cell2struct(reshape(getappdata(0,ad.tmp_storage), ad.size), ad.fields, 1);
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
  % (c) E.Farhi, ILL. License: EUPL.
  
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
      'Tag', options.Tag, ...
      'CloseRequestFcn', ...
        [ 'setappdata(0,''' tmp_storage ''', get(get(gcbf,''Children''),''Data'')); delete(gcbf)' ]);
  
  % override default mechanism for Table update when closing
  if isfield(options, 'CloseRequestFcn') && ~strcmp(options.CreateMode, 'modal')
    set(f, 'CloseRequestFcn', options.CloseRequestFcn);
  end
  
  % create the Table content to display, handle array of structures
  % check if the structure values are numeric, logical, or char
  Data0       = cell(numel(fields), numel(structure));
  fields_type = cell(size(Data0));
  for index_s = 1:numel(structure)
    Data0(:,index_s) = struct2cell(structure(index_s));
    for index=1:numel(fields)
      item = Data0{index,index_s};
      if ~isnumeric(item) && ~islogical(item) && ~ischar(item)
        fields_type{index,index_s} = class(item);
        Data0{index,indep_s} = class2str('', Data0{index,index_s}, 'eval');
      end
    end
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
  
  % determine the window size to show
  TextWidth = 12;
  TextHeight= 30;
  % height is given by the number of fields
  height = (numel(options.ListString)+3)*TextHeight;
  % width is given by the length of the longest RowName
  width = max(cellfun(@numel,options.ListString))*TextWidth + 10*options.FontSize;
  % compare to current window size
  p = get(f, 'Position');
  p(3) = width;
  if p(4) > height, p(4) = height; end
  set(f, 'Position',p);

  % create the table
  t = uitable('Parent',f, ...
    'RowName',options.ListString, 'Data', Data0, ...
    'ColumnEditable',true, ...
    'FontSize',options.FontSize, ...
    'Units','normalized', 'Position', [0.05 0.2 .9 .7 ], ...
    'ColumnWidth','auto','TooltipString',options.TooltipString);
    
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
    
    if numel(Data) ~= numel(Data0), return; end
    % restore initial structure field type when set
    Data = reshape(Data, ad.size);
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

function structure = structdlg(structure,options)
  % allow to modify a structure in a table
  %
  % input:
  %   structure: the initial struct to edit
  %   options:   a set of options, namely:
  %     options.Name: the name of the dialogue
  %     options.FontSize: the FontSize used to display the table
  %       default: use the system FontSize.
  %     options.ListString: the labels to be used for each structure field. 
  %       default: use the structure member names
  %     options.TooltipString: a string to display (as help)
  
  Data0   = struct2cell(structure);   % initial Data
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
    options.ListString = fields;
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
  setappdata(0,tmp_storage,[]);   % we shall collect the Table content here
  
  % create a uitable from the structure fields and values
  f = figure('Name',options.Name, 'MenuBar','none', ...
      'CloseRequestFcn', ...
        [ 'setappdata(0,''' tmp_storage ''', get(get(gcbf,''Children''),''Data'')); delete(gcbf)' ]);
        
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

  t = uitable('Parent',f, ...
    'RowName',options.ListString, 'Data', Data0, ...
    'ColumnEditable',true, ...
    'FontSize',options.FontSize, ...
    'Units','normalized', 'Position', [0.05 0.2 .9 .7 ], ...
    'ColumnWidth','auto','TooltipString',options.TooltipString);
  
  % when the figure is deleted, we should get the uitable Data back
  % wait for figure close
  uiwait(f);
  Data = getappdata(0,tmp_storage);
  if numel(Data) ~= numel(Data0), return; end
  % assemble new options...
  
  structure = cell2struct(Data, fields, 1);
  rmappdata(0, tmp_storage);

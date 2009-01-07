function [hIcon,config,object]=iView_private_icon(instance, action, object, object_index)
% iView_private_icon_new Add a new data set inside the iView instance
% The object may be an iData, or any other type, which is then sent to iData
% before importation.
% action: 
% load, check, resize, documents, open, select_all, deselect_all, delete, saveas
% check_load

hIcon=[];

if isempty(action), action='check'; end

% actions which do not require an object to be given
switch action
  case 'resize'     % re-arrange icons in instance window
    if isempty(object), object = getappdata(instance, 'Data'); end
    iView_private_icon(instance, 'check', object);
    return
  case 'documents'  % populate the Documents menu
    iView_private_icon_documents(instance);
    return
  case {'select_all','deselect_all'}
    objects = get(gcf,'Children');
    for index=1:length(objects)
      this=objects(index);
      if strcmp(get(this,'type'),'uicontrol')
        if ~isempty(strmatch(get(this,'Style'), {'togglebutton','checkbox','pushbutton'}))
          set(this, 'Value', object);
        end
      end
    end
    return
end % switch action

if nargin < 3, return; end

Data = getappdata(instance, 'Data');
if isempty(object), return; end
if ischar(object) % special data selections
  switch object
  case 'selection'
    object = iView_private_selection(instance);
  case 'all'
    object = Data;
  end
end

config = iView_private_config(instance, 'load');
% make sure we have an iData
if ~isempty(object) & ischar(object)
  object=cellstr(object);
end
if ~isempty(object) & iscellstr(object)
  new_object=[];
  for index=1:length(object)
    this = findobj(Data, 'Tag', object{index});
    if isempty(this)
      % object has no corresponding iData
      % remove uicontrol and exit
      hIcon = findobj(instance, 'Tag', object{index});
      delete(hIcon);
    else
      new_object= [new_object this ];
    end
  end
  object=new_object;
  if isempty(object), return; end
end
if ~isa(object, 'iData')
  object=iData(object); % will also popup FileSelector from iLoad/uigetfiles if object=''
end

% actions which handle iData array
switch action
case 'open'
  hIcon=subplot(object);
  return
case 'save'
  if ~isempty(object)
    format=config.OutputFormat;
    [hIcon,format]=saveas(object, format);    
    if ~isempty(hIcon)
      disp('Objects')
      disp(object)
      disp('into')
      disp(hIcon);
    end
  end
  return
case 'saveas'
  if ~isempty(object)
    [hIcon,format]=saveas(object, 'gui');
    config.OutputFormat = format;
    config = iView_private_config(instance, 'save');
    
    if ~isempty(hIcon)
      disp('Objects')
      disp(object)
      disp('into')
      disp(hIcon);
    end
  end
  return
case 'delete'
  tokeep = 1:length(Data);
  DataTags = get(Data, 'Tag');
  for index=1:length(object)
    this = ind2sub(object, index);
    this_delete = strmatch(this.Tag, DataTags, 'exact');
    hIcon = findobj(instance,'Type','uicontrol','Style', config.IconStyle, 'Tag', this.Tag);
    delete(hIcon);
    tokeep(this_delete) = 0;
  end
  Data = ind2sub(Data, tokeep(find(tokeep)));
  setappdata(instance, 'Data', Data);
  iView_private_icon(instance, 'documents', []);
  iView_private_icon(instance, 'check', []);
  return
end

% handle iData arrays
if length(object) > 1
  for index=1:length(object(:))
    [i,c,d]= iView_private_icon(instance, action, object(index), index);
    object{index}=d;
    hIcon= [hIcon i ];
  end
  return
end

% actions which handle single iData 
if isempty(object), return; end

index  = strmatch(object.Tag, get(Data, 'Tag'), 'exact');
hIcon  = findobj(instance, 'Tag', object.Tag);

switch action
case {'check','check_load'}
  % remove uicontrols that do not have corresponding iData
  hIcon = findobj(instance,'Type','uimenu', 'Tag', [ 'Doc_' object.Tag ]);
  if ~isempty(hIcon) & isempty(index)
    delete(hIcon); hIcon=[];
  end
  hIcon = findobj(instance,'Type','uicontrol','Style', config.IconStyle, 'Tag', object.Tag);
  if length(hIcon) ~= length(index)
    delete(hIcon); hIcon=[];
  end
  if nargin > 3 && isempty(index), index=object_index; end
case 'load'
  % check if the element already exists, add it if no
  if ~isempty(index) % already exists 
    % make a copy of the object
    object = copyobj(object);
  end
  Data = [ Data object ]; % store new object
  setappdata(instance, 'Data', Data);
  index= length(Data);
case 'properties'
  a = object;
  T=a.Title; if iscell(T), T=T{1}; end
  T=strtrim(T);
  cmd = char(a.Command{end});
  properties={ [ 'Data ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ], ...
               [ 'Title: "' T '"' ], ...
               [ 'Source: ' a.Source ], ...
               [ 'Date: ' a.Date ], ...
               [ 'Last command: ' cmd ]};
  if ~isempty(a.Label)
    properties{end+1} = [ 'Label: ' a.Label ];
  end
  if length(getaxis(a))
    properties{end+1} = '[Rank]         [Value] [Description]';
    for index=0:length(getaxis(a))
      [v, l] = getaxis(a, num2str(index));
      x      = getaxis(a, index);
      m      = get(a, 'Monitor');
      if index==0 & not(all(m==1) | all(m==0))
        m      = get(a, 'Monitor');
        properties{end+1} = sprintf('%6i %15s  %s [%g:%g] (per monitor)', index, v, l, min(x(:)./m(:)), max(x(:)./m(:)));
      else
        properties{end+1} = sprintf('%6i %15s  %s [%g:%g]', index, v, l, min(x(:)), max(x(:)));
      end
    end
  end
  S=a.Source;
  [dummy, fS] = fileparts(S);
  titl ={ T ; [ a.Tag ' <' fS '>' ]};
  if length(T) > 23, T=[ T(1:20) '...' ]; end
  if length(S) > 23, S=[ '...' S(end-20:end) ]; end
  if length(cmd) > 23, cmd = [ cmd(1:20) '...' ]; end

  cdata = get(hIcon, 'cdata');
  if ~isempty(cdata)
    msgbox(properties, ['iView: Properties of ' T ' <' S '>'],'custom', cdata);
  else
    msgbox(properties, ['iView: Properties of ' T ' <' S '>']);
  end
  return
otherwise
end

if isempty(index), return; end

% determine position of icon: we use the icon size from the Config
% position is set from iData index (end of array), icon size and window width
[iRow,iColumn]  = iView_private_icon_position(instance, index);

NL = sprintf('\n');
tooltip = [ object.Title NL object.Source NL 'Tag ' object.Tag '; Data size [' num2str(size(object)) ']' NL ];
for ax_index = 0:ndims(object)
  [axisdef, lab] = getaxis(object, num2str(ax_index));
  lab=strtrim(lab);
  if ax_index==0, 
    if isempty(lab), tooltip = [ tooltip 'Signal' ]; else tooltip = [ tooltip lab ]; end
  elseif ~isempty(lab)
    tooltip = [ tooltip ' vs ' lab ]; 
  end
end

if strcmp(action, 'load') | strcmp(action, 'check_load')
  if isempty(hIcon)
    % create icon there
    % we can not use the Callback for the button, else only drag n drop feature can not be handled (inactive)
    % for this we set the Enable to inactive. This unactivates Callback, but we can manage
    % all events, e.g. ButtonDownFcn
    % The Enable=inactive removes the Tooltip
    if ~isfield(object.UserData, 'cdata')
      frame = getframe(object, config.IconSize);
      cdata = frame.cdata;
      object.UserData.cdata=cdata;
    else
      cdata = object.UserData.cdata;
    end
    hIcon = uicontrol(instance,'Style', config.IconStyle,...
                    'Tag', object.Tag, ...
                    'Value',0, ...
                    'Enable', 'on',...
                    'CData', cdata);
    if ~strcmp(config.DragDropNoTooltips,'no') & config.DragDropNoTooltips
      set(hIcon, 'Enable','inactive', 'ButtonDownFcn', 'iview(gcf, ''mouse_down'', gco);');
    end
  end
else
  hIcon = findobj(instance, 'Tag', object.Tag);
end
if isempty(hIcon), return; end

% icon Tag must be the one from the iData object
set(hIcon, ... 
  'String',object.Title,...
  'Position',[iColumn iRow 2*config.IconSize config.IconSize], ...
  'ToolTipString', tooltip, ...
  'Tag', object.Tag);

% uicontext menu on icons
%   select
%   deselect
%   open
%   save
%   rename
%   cut
%   copy
%   paste (merge)
%   delete (close)
%   properties
%   edit (file)
%   edit (axes, signal, title, ...)
figure(instance);
cmenu = uicontextmenu('Parent',instance);
set(hIcon, 'UIContextMenu', cmenu);
uimenu(cmenu, 'Label', object.Title);
uimenu(cmenu, 'Label', 'Open', 'Separator','on', 'Callback', [ 'set(gco,''Value'',1); iview(gcf, ''open_data'', ''selection'');' ]);
uimenu(cmenu, 'Label', 'Save', 'Callback', [ 'set(gco,''Value'',1); iview(gcf, ''save_data'', ''selection'');' ]);
uimenu(cmenu, 'Label', 'Save as...', 'Callback', [ 'set(gco,''Value'',1); iview(gcf, ''saveas_data'', ''selection'');' ]);
uimenu(cmenu, 'Label', 'Close', 'Callback', [ 'set(gco,''Value'',1); iview(gcf, ''close_data'', ''selection'');' ]);
uimenu(cmenu, 'Label', 'Properties', 'Callback', [ 'set(gco,''Value'',1); iview(gcf, ''properties'', ''selection'');' ], 'Separator','on');

% =============================================================================
%                    INLINE functions
% =============================================================================
% iView_private_icon_position
% iView_private_icon_documents

function [iRow,iColumn]  = iView_private_icon_position(instance, index)
% iView_private_icon_position determine position of the indexed element in an 
% iView instance. 

  hIcon = [];

  config = iView_private_config(instance, 'load');

  % size of Icon is [ config.IconWidth*config.IconSize config.IconSize ]

  % number of elements per row
  instancePosition = get(instance, 'Position');
  instanceWidth    = instancePosition(3);
  nbElementsPerRow = max(1,floor(instanceWidth/((config.IconWidth+.25)*config.IconSize)));

  if instanceWidth < (config.IconWidth+.25)*config.IconSize
    instancePosition(3) = (config.IconWidth+.25)*config.IconSize;
    set(instance, 'Position', instancePosition);
    instanceWidth    = instancePosition(3);
    figureslider(instance, 'resize');
  end

  % get icon index
  iColumn = mod(index-1, nbElementsPerRow);
  iRow    = floor((index-1)/nbElementsPerRow);
  % make it pixels in the window (with a gap of 1/4 between elements, and 1/2 shift)
  iRow    = iRow    *(1.25*config.IconSize) + config.IconSize/4;
  iColumn = iColumn *((config.IconWidth+.25)*config.IconSize) + config.IconSize/4;
  % fill in from the top (revert iRow)
  iRow = instancePosition(4) - iRow - config.IconSize;



function iView_private_icon_documents(instance)
% iView_private_icon_documents Update Documents menu with list of loaded iData objects

  documents = findobj(instance, 'Tag', 'Documents');
  items     = get(documents, 'Children');
  
  % clear Static handles from the list (not from menu)
  index     = strmatch('Static', get(items, 'Tag'), 'exact');
  items(index)=[];
  % delete other entries (as we shall rebuild menu items)
  delete(items);
  
  % populate Documents
  Data = getappdata(instance, 'Data');
  
  for index=1:length(Data)
    if length(Data) == 1, this=Data; else this=Data(index); end
    h = uimenu(documents, 'Label', Data(index).Title, ...
        'Callback', [ 'iview(gcf, ''open_data'', ''' this.Tag ''');' ], ...
        'Tag',[ 'Doc_'  this.Tag ]);
    if index==1
      set(h, 'Separator', 'on');
    end
  end
  

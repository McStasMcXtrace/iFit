function hIcon=iView_private_icon(instance, action, object)
% iView_private_icon_new Add a new data set inside the iView instance
% The object may be an iData, or any other type, which is then sent to iData
% before importation.
% action: load, check, resize, documents, open

hIcon=[];
if nargin < 3, return; end

if isempty(action), action='check'; end

if strcmp(action, 'resize')     % re-arrange icons in instance window
  iView_private_icon_resize(instance);
  return
end
if strcmp(action, 'documents')  % populate the Documents menu
  iView_private_icon_documents(instance);
  return
end

if isempty(object), return; end

Data = getappdata(instance, 'Data');
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

% handle iData arrays
if length(object) > 1
  for index=1:length(object(:))
    hIcon= [hIcon iView_private_icon(instance, action, object(index))];
  end
  return
end

config = iView_private_config(instance, 'load');
index  = strmatch(object.Tag, get(Data, 'Tag'), 'exact');

switch action
case 'check'
  hIcon = findobj(instance, 'Tag', object.Tag);
  if ~isempty(hIcon) & isempty(index)
    delete(hIcon);
  end
case 'open'
  hIcon=subplot(object);
  return
case 'delete'
  if ~isempty(index)
    % remove iData and corresponding uicontrol
    if length(Data) > 1, 
      Data(index) = [];
      Data=Data(find(~isempty(Data)));
    else Data=[]; end
    setappdata(instance, 'Data', Data);
    iView_private_icon(instance, 'check', object);
  end
case 'saveas'
  if ~isempty(object)
    hIcon=saveas(object, 'gui');
    if ~isempty(hIcon)
      fprintf(1, 'Object %s\n  saved into file %s\n', char(object), hIcon);
    end
  end
  return
case 'load'
  % check if the element already exists, add it if no
  if ~isempty(index) % already exists 
    % make a copy of the object
    object = copyobj(object);
  end
  Data = [ Data object ]; % store new object
  setappdata(instance, 'Data', Data);
  index= length(Data);
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
  lab=deblank(lab);
  if ax_index==0, 
    if isempty(lab), tooltip = [ tooltip 'Signal' ]; else tooltip = [ tooltip lab ]; end
  elseif ~isempty(lab)
    tooltip = [ tooltip ' vs ' lab ]; 
  end
end

if strcmp(action, 'load')
  % create icon there
  frame = getframe(object, config.IconSize);
  hIcon = uicontrol(instance,'Style','togglebutton',...
                  'Tag', object.Tag, ...
                  'Value',0, ...
                  'Callback', 'iview(gcf, ''mouse_down'', gco);', ...
                  'CData', frame.cdata);
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

% install mouse event handling callbacks (on icons only !!!)
%   move icon (re-order)
%   drag-drop
%   double click

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
uimenu(cmenu, 'Label', 'Open', 'Separator','on', 'Callback', [ 'iview(gcf, ''open_data'', ''' object.Tag ''');' ]);
uimenu(cmenu, 'Label', 'Save as...', 'Callback', [ 'iview(gcf, ''save_data'', ''' object.Tag ''');' ]);
uimenu(cmenu, 'Label', 'Close', 'Callback', [ 'iview(gcf, ''close_data'', ''' object.Tag ''');' ]);

% =============================================================================
%                    INLINE functions
% =============================================================================
% iView_private_icon_position
% iView_private_icon_resize
% iView_private_icon_documents

function [iRow,iColumn]  = iView_private_icon_position(instance, index)
% iView_private_icon_position determine position of the indexed element in an 
% iView instance. 

  hIcon = [];

  config = iView_private_config(instance, 'load');

  % size of Icon is [ 2*config.IconSize config.IconSize ]

  % number of elements per row
  instancePosition = get(instance, 'Position');
  instanceWidth    = instancePosition(3);
  nbElementsPerRow = max(1,floor(instanceWidth/(2.25*config.IconSize)));

  if instanceWidth < 2.25*config.IconSize
    instancePosition(3) = 2.25*config.IconSize;
    set(instance, 'Position', instancePosition);
    instanceWidth    = instancePosition(3);
    figureslider(instance, 'resize');
  end

  % get icon index
  iColumn = mod(index-1, nbElementsPerRow);
  iRow    = floor((index-1)/nbElementsPerRow);
  % make it pixels in the window (with a gap of 1/4 between elements, and 1/2 shift)
  iRow    = iRow    *(1.25*config.IconSize) + config.IconSize/4;
  iColumn = iColumn *(2.25*config.IconSize) + config.IconSize/4;
  % fill in from the top (revert iRow)
  iRow = instancePosition(4) - iRow - config.IconSize;

function iView_private_icon_resize(instance)
% iView_private_icon_resize reposition icons inside iView instance

  Data = getappdata(instance, 'Data');
  iView_private_icon(instance, 'check', Data);

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
    if index==1
      uimenu(documents, 'Label', Data(index).Title, 'Callback', [ 'iview(gcf, ''open_data'', ''' this.Tag ''');' ], 'Separator', 'on');
    else
      uimenu(documents, 'Label', Data(index).Title, 'Callback', [ 'iview(gcf, ''open_data'', ''' this.Tag ''');' ]);
    end
  end
  

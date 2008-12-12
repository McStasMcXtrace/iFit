function iView_private_icon(instance, action, object)
% iView_private_icon_new Add a new data set inside the iView instance
% The object may be an iData, or any other type, which is then sent to iData
% before importation

if nargin < 3, return; end

if isempty(action), action='check'; end

if strcmp(action, 'resize')
  iView_private_icon_resize(instance);
  return
end

if isempty(object), return; end

% make sure we have an iData
if ~isempty(object) & ~isa(object, 'iData')
  object=iData(object);
end

% handle iData arrays
if length(object) > 1
  for index=1:length(object(:))
    iView_private_icon(instance, action, object(index));
  end
  return
end

% check if the element already exists, add it if no
Data = getappdata(instance, 'Data');
config = getappdata(instance, 'Config');
index = get(Data, 'Tag');

index = strmatch(object.Tag, index, 'exact');

if isempty(index) | strcmp(action, 'load')
  if ~isempty(index)
    % make a copy of the object
    object = copyobj(object);
  end
  Data = [ Data object ]; % store new object
  setappdata(instance, 'Data', Data);
  index= length(Data);
end

% determine position of icon: we use the icon size from the Config
% position is set from iData index (end of array), icon size and window width
[iRow,iColumn]  = iView_private_icon_position(instance, index);

NL = sprintf('\n');
tooltip = [ object.Title NL object.Source NL 'Data size [' num2str(size(object)) ']' NL ];
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

% =============================================================================
%                    INLINE functions
% =============================================================================
% iView_private_icon_position

function [iRow,iColumn]  = iView_private_icon_position(instance, index)
% iView_private_icon_position determine position of the indexed element in an 
% iView instance. 

hIcon = [];

config = getappdata(instance, 'Config');

% size of Icon is [ 2*config.IconSize config.IconSize ]

% number of elements per row
instancePosition = get(instance, 'Position');
instanceWidth    = instancePosition(3);
nbElementsPerRow = ceil(instanceWidth/(3*config.IconSize));

if instanceWidth < 2.25*config.IconSize
  instancePosition(3) = 2.25*config.IconSize;
  set(instance, 'Position', instancePosition);
  instanceWidth    = instancePosition(3);
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
% reposition icons inside iView instance

Data = getappdata(instance, 'Data');
iView_private_icon(instance, 'check', Data);



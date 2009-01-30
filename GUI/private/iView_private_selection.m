function [selection, selectedIndex, selected] = iView_private_selection(instance, tags)
% get selected Data sets (all Value=1 objects). selected=handle of icons
% returns: 
% selection    = iData array (only selected elements), 
% selectedIndex= index of these elements in the window iData array
% selected     = handle to the selected icons (uicontrols)

config = iView_private_config(instance, 'load');
Data        = getappdata(instance, 'Data');
DataTags    = get(Data, 'Tag');     % all Data tags in window
selection   = []; selectedIndex = []; selected = [];

if nargin < 2
  tags    = findobj(instance,'Type','uicontrol','Style', config.IconStyle, 'Value',1);
  tags    = get(tags,'Tag');
elseif isa(tags,'iData') | ishandle(tags)
	tags    = get(tags,'Tag');
end
if isempty(tags), return; end
if ~iscellstr(tags), tags=cellstr(tags); end
selectedIndex=zeros(size(tags));
	
% get index of these selected icons in Data
for index=1:length(tags)
  i = strmatch(tags{index}, DataTags, 'exact');
  if length(i) > 1
    warning([ 'Duplicated icon ' tags{index} ' in figure ' num2str(instance) ])
  end
  if ~isempty(i), 
    selectedIndex(index)=i(end);
  end
end
selectedIndex = selectedIndex(find(selectedIndex));

selection     = ind2sub(Data, selectedIndex);              % selected iData objects (initial)
  


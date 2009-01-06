function [selection, selectedIndex, selected] = iView_private_selection(instance)
% get selected Data sets (all Value=1 objects). selected=handle of icons
% returns: 
% selection    = iData array (only selected elements), 
% selectedIndex= index of these elements in the window iData array
% selected     = handle to the selected icons (uicontrols)
  config = iView_private_config(instance, 'load');
  selected    = findobj(instance,'Type','uicontrol','Style', config.IconStyle, 'Value',1);
  Data        = getappdata(instance, 'Data');
  DataTags    = get(Data, 'Tag');     % all Data tags in window
  selectedIndex=zeros(size(selected));
  % get index of these selected icons in Data
  for index=1:length(selected)
    i = strmatch(get(selected(index), 'Tag'), DataTags, 'exact');
    if length(i) > 1
      warning([ 'Duplicated icon ' get(selected(index), 'Tag') ' in figure ' num2str(instance) ])
    end
    if ~isempty(i), 
      selectedIndex(index)=i(end);
    end
  end
  selectedIndex = selectedIndex(find(selectedIndex));
  
  selection     = iView_private_icon_idata(Data, selectedIndex);              % selected iData objects (initial)
    


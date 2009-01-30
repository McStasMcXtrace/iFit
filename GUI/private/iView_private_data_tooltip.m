function tooltip = iView_private_data_tooltip(object)
% iView_private_data_tooltip creates a tool tip for an icon, based on data
NL = sprintf('\n');
tooltip = object.Title;
if ~isempty(object.Label)
	tooltip = [ tooltip ' (' object.Label ')' ];
end
tooltip = [ tooltip NL object.Source NL ...
	 'Tag ' object.Tag '; Data size [' num2str(size(object)) ']' NL ];
for ax_index = 0:ndims(object)
  [axisdef, lab] = getaxis(object, num2str(ax_index));
  lab=strtrim(lab);
  if ax_index==0, 
    if isempty(lab), tooltip = [ tooltip 'Signal' ]; else tooltip = [ tooltip lab ]; end
  elseif ~isempty(lab)
    tooltip = [ tooltip ' vs ' lab ]; 
  end
end

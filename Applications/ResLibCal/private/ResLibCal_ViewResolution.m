function [out, h] = ResLibCal_ViewResolution(out, dim)
% ResLibCal_ViewResolution: open the Resolution 2D/3D plot view
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  h = findobj(0, 'Tag',[ 'ResLibCal_View' num2str(dim)]);
  if isempty(h)
    if dim~=1, name=sprintf('(%iD)', dim); else name='Matrix'; end
    h = figure('Name',[ 'ResLibCal: View Resolution ' name ], ...
               'Tag', [ 'ResLibCal_View' num2str(dim)], 'ToolBar','figure');
    p = get(h, 'Position'); p(3:4) = [ 640 480 ]; set(h, 'Position',p);
  else
    figure(h);
  end

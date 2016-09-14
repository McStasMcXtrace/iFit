% ResLibCal_fig: finds current ResLibCal main GUI, as well as other
%   controls.
%   A cache is built when accessing the controls for the first time, and
%   is re-used afterwards for faster access (e.g. Matlab >= R2014b)
%
%   ResLibCal_fig:      returns the main GUI handle
%   ResLibCal_fig(tag): returns the 'Tag' specified
function f=ResLibCal_fig(tag)

  persistent fig handles

  if ~ishandle(fig), fig=[]; end
  if isempty(fig)
    fig = findall(0, 'Tag','ResLibCal');
    if length(fig) > 1, delete(fig(2:end)); end
    handles = [];
  end
  
  if isempty(handles), handles.fig = fig; end

  if nargin == 0
    f=fig;
  else
    if ~isfield(handles, tag) || ~ishandle(handles.(tag))
      handles.(tag) = findobj(fig, 'Tag', tag);
    end
    f = handles.(tag);
  end

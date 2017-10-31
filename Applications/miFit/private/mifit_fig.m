function f=mifit_fig(tag)
% [internal] mifit_fig: search for a given Tag in Application or main Figure if ommitted.
  persistent fig
  persistent handles
  
  if ~ishandle(fig), fig=[]; end
  if isempty(fig)
    fig = findall(0, 'Tag','miFit');
    if length(fig) > 1, delete(fig(2:end)); end % unique instance
    handles = [];
  end
  if isempty(handles)
    handles.fig = fig;
  end

  if nargin == 0
    f=fig;
  elseif ~isempty(handles)
    if strcmp(tag,'handles'), f=handles; return;
    elseif ~isfield(handles, tag) || (~isempty(handles.(tag)) && any(~ishandle(handles.(tag)) ))
      handles.(tag) = []; end
    if isempty(handles.(tag))
      handles.(tag) = findall(fig, 'Tag', tag);
      if isempty(handles.(tag))
        handles.(tag) = findall(0, 'Tag', tag);
      end
    end
    f = handles.(tag);
  else f = [];
  end

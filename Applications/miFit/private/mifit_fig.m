function f=mifit_fig(tag)
% search for a given Tag in Application or main Figure if ommitted.
  persistent fig handles
  
  if ~ishandle(fig), fig=[]; end
  if isempty(fig)
    fig = findall(0, 'Tag','miFit');
    if length(fig) > 1, delete(fig(2:end)); end % unique instance
    handles = [];
  end

  if nargin == 0
    f=fig;
  else
    if ~isfield(handles, tag) || ~ishandle(handles.(tag)) 
      handles.(tag) = []; end
    if isempty(handles.(tag))
      handles.(tag) = findobj(fig, 'Tag', tag);
      if isempty(handles.(tag))
        handles.(tag) = findall(0,'Tag', tag);
      end
    end
    f = handles.(tag);
  end

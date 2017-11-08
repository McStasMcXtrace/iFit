function selected = mifit_Optimizers_Set(varargin)
% [internal] mifit_Optimizers_Set: set current optimizer, and check its label in the menu
  % set the optimizer to use
  selected0 = getappdata(mifit_fig, 'CurrentOptimizer');
  selected  = selected0;
  if nargin && ~isempty(varargin{1})
    if ishandle(varargin{1})
      selected = get(varargin{1},'UserData');
    elseif ischar(varargin{1})
      selected = varargin{1};
    end
  end
  if isempty(selected), selected='fmin'; end  % default when not set
  if ~strcmp(selected, selected0) % when we change the optimizer, we reset its Configuration
    options = feval(selected, 'defaults');
    setappdata(mifit_fig,'CurrentOptimizerConfig', options);
  end
  mifit_disp([ '[Optimizer] Setting optimizer to "' selected '"' ]);
  setappdata(mifit_fig,'CurrentOptimizer', selected);

  % we get the CurrentOptimizer, and check it. uncheck the others
  hmodels = mifit_fig('Menu_Optimizers');
  hoptims = get(hmodels,'Children');
  names   = get(hoptims,'UserData');
  index   = strcmp(selected, names);
  set(hoptims,'Checked','off');
  set(hoptims(index),'Checked','on');

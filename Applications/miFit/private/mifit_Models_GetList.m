function [ifuncs, labels,indices,handles] = mifit_Models_GetList(varargin)
% [internal] mifit_Models_GetList: get the list of iFunc models (not expressions) in the Models menu
% indices is the index of 'static' Models in the whole list.
  models = getappdata(mifit_fig,'Models');
  ifuncs = []; labels = {}; indices = []; handles = [];
  for index=1:numel(models)
    this = models{index};
    if ~isempty(this.callback) && isa(this.callback, 'iFunc')
      ifuncs = [ ifuncs this.callback ];
      handles= [ handles this.handle ];
      labels{end+1} = this.label;
      indices(end+1) = index;
    end
  end

function mifit_List_Model_push(d, flag_replace)
% [internal] mifit_List_Model_push: put/replace a new model in the Models menu
  if isempty(d),       return; end
  if nargin == 1, flag_replace = []; end
  if iscell(d)
    for index=1:numel(d)
      mifit_List_Model_push(d{index}, flag_replace);
    end
    return
  end
  if isa(d, 'iData')
    mifit_List_Data_push(d, flag_replace);
    return
  end
  fig = mifit_fig;

  % update AppData Stack
  if numel(d) > 1, d = d(:); end
  Models = getappdata(fig, 'Models');

  if strcmp(flag_replace,'replace')
    % we search for Models that have the same Tag, and replace them
    for index=1:numel(d)
      for mindex = 1:numel(Models)
        obj = Models{mindex};
        if numel(d) == 1, this = d; else this = d(index); end
        if isa(obj.object, 'iFunc') && strcmp(this.Tag, obj.object.Tag)
          obj.object = this;
          Models{mindex} = obj;
        end
      end
    end
    setappdata(fig, 'Models', Models);
  else
    mifit_Models_Add_Entry(d);
    return
  end

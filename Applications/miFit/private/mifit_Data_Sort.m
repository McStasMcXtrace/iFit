function mifit_Data_Sort(varargin)
% Data/Sort: sort selected data sets or all when none selected
%   sort can be: ID/Tag, file, name, creation date, intensity, ...

  % get the selected data or the whole list
  [d, index_selected] = mifit_List_Data_pull;
  D=getappdata(mifit_fig, 'Data');
  if isempty(index_selected) || numel(d) == 0, 
    d=D;
    index_selected=1:numel(D);
  end

  % do not sort when nothing there
  if isempty(index_selected) || numel(d) ==0, return; end
  
  % display a dialogue to select sorting criteria
  op = {'iD number','Data size', 'Intensity (sum)', 'Title/Name/Label/User','File', ...
    'Date', ...
    'Move data sets at the beginning', ...
    'Move data sets at the end', ...
    'Invert selection (select those unselected)', ...
    'Revert order (first goes to last)'};
  PromptString = cellstr(char([ d{:} ]));
  if numel(PromptString) > 3
    PromptString=[ PromptString(1:3) ; ' ...' ];
  end
  [selection, ok] = listdlg('ListString', op, 'ListSize',[500 160], ...
    'Name',[ 'miFit: Select sorting criteria to apply on ' ...
      num2str(numel(index_selected)) ' data set(s)' ], ...
      'PromptString', PromptString, 'SelectionMode','single');
      
  if isempty(selection) || ok ~= 1, return; end
  index = [];
  selection = op{selection};
  
  % Sorting than can apply on only one data set:
  %   'Move data sets at the beginning'
  %   'Move data sets at the end'
  %   'Invert selection (select those unselected)'
  % Other actions require more than 1 data set to be selected. We thus extend to
  % the whole list when only 1 selected.
  if ~any(strcmp(strtok(selection), { 'Move','Invert'})) && numel(d) < 2
    d=D;
    index_selected=1:numel(D);
  end
  
  if iscell(d), d=[ d{:} ]; end
  
  % sort
  switch selection
  case 'iD number'
    criteria = get(d,'Tag');
  case 'Data size'
    criteria = zeros(size(d));
    for index=1:numel(d)
      criteria(index) = prod(size(d(index)));
    end
    index=[];
  case 'Intensity (sum)'
    criteria = sum(d,0);
  case 'Title/Name/Label/User'
    criteria = strtrim(strcat(...
      get(d,'Title')', ' ' , ...
      get(d,'Name')', ' ' , ...
      get(d,'Label')', ' ' , ...
      get(d,'User')'));
  case 'File'
    criteria = get(d,'Source');
  case 'Date'
    criteria = get(d,'Date');
  case 'Move data sets at the end'
    selected      = d; % == D(index_selected)
    index         = ones(1,numel(D));
    index(index_selected) = 0;
    index         = find(index);
    not_selected  = D(index);
    D = { not_selected{:} selected{:} };
    setappdata(mifit_fig, 'Data',D);
    mifit_List_Data_UpdateStrings;
    mifit_Edit_Select_All([], (numel(not_selected)+1):numel(D));
    mifit_History_push;
    return
  case 'Move data sets at the beginning'
    selected      = d; % == D(index_selected)
    index         = ones(1,numel(D));
    index(index_selected) = 0;
    index         = find(index);
    not_selected  = D(index);
    D = { selected{:} not_selected{:} };
    setappdata(mifit_fig, 'Data',D);
    mifit_List_Data_UpdateStrings;
    mifit_Edit_Select_All([], 1:numel(selected));
    mifit_History_push;
    return
  case 'Invert selection (select those unselected)'
    index = ones(1,numel(D));
    index(index_selected) = 0;
    index = find(index);
    mifit_Edit_Select_All([], index);
    return
  case 'Revert order (first goes to last)'
    criteria = numel(d):-1:1;
  otherwise
    disp([ mfilename ': unsupported sorting criteria ' selection ]);
    return
  end
  
  % update the Data set List
  if isempty(index)
    if ~iscell(criteria) && ~isnumeric(criteria)
      return;
    end
    [criteria,index] = sort(criteria);
  end

  D(index_selected) = D(index_selected(index)); % re-order
  setappdata(mifit_fig, 'Data',D);
  mifit_List_Data_UpdateStrings;
  mifit_Edit_Select_All([], sort(index_selected(index)));
  mifit_History_push;
  
  % reset the CurrentDataSet
  
  
  
  

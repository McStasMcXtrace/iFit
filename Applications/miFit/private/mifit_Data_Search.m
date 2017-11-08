function mifit_Data_Search(varargin)
% Data/Search: search all Data sets for a token
  prompt = {'Enter the item to search for (string, case sensitive):'};
  name   = 'miFit: Search in all data sets';
  numlines=1;
  answer =inputdlg(prompt,name,numlines);
  if isempty(answer),    return; end  % cancel
  if isempty(answer{1}), return; end  % empty pattern
  
  D=getappdata(mifit_fig, 'Data'); D = [ D{:} ];
  
  mifit_disp([ '[Data search] Data sets mentioning "' answer{1} '" ...' ]);
  set(mifit_fig,'Pointer','watch');
  % search for token in strings
  match = strfind(D, answer{1}, 'case'); match=match(:);
  index_selected = ~cellfun(@isempty, match);
  % search for token in field names
  match = findfield(D, answer{1}, 'case'); match=match(:);
  match = ~cellfun(@isempty, match);
  index_selected = index_selected | match;
  % search for tokens in the char(iData) string (as displayed)
  match          = strfind(cellstr(char(D)),answer{1}); match=match(:);
  match = ~cellfun(@isempty, match);
  index_selected = index_selected | match;
  index_selected = find(index_selected);
  if isempty(index_selected)
    mifit_disp('None');
  else
    mifit_disp(char(D(index_selected)));
    mifit_Edit_Select_All([], index_selected);
  end
  set(mifit_fig,'Pointer','arrow');

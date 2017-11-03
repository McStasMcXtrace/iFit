function mifit_List_Data_UpdateStrings
% [internal] mifit_List_Data_UpdateStrings: update the List labels
  fig = mifit_fig;
  
  Data = getappdata(fig, 'Data');
  
  hObject        = mifit_fig('List_Data_Files');
  list           = {};
  if numel(Data) > 1
    for index=1:numel(Data)
      [~,label] = strtok(char(Data{index}));
      % we remove the initial 'iData' word
      list{end+1} = strtrim(label);
    end
  else
    list{end+1} = char(Data{1});
  end
  set(hObject,'String', list, 'Value', []);

function iView_private_documents(instance)
% iView_private_documents Update Update Documents menu with list of loaded iData objects
% and list of labels

  documents = findobj(instance, 'Tag', 'Documents');
  items     = get(documents, 'Children');
  
  % clear Static handles from the list (not from menu)
  index     = strmatch('Static', get(items, 'Tag'), 'exact');
  items(index)=[];
  % delete other entries (as we shall rebuild menu items)
  delete(items);
  
  % populate Documents
  Data = getappdata(instance, 'Data');
  
  for index=1:length(Data)
    if length(Data) == 1, this=Data; else this=Data(index); end
    h = uimenu(documents, 'Label', this.Title, ...
        'Callback', [ 'iview(gcf, ''data_open'', ''' this.Tag ''');' ], ...
        'Tag',[ 'Doc_'  this.Tag ]);
    if index==1
      set(h, 'Separator', 'on');
    end
  end
  
  % update list of known labels in the 'Set Label as' item
  setlabelas = findobj(documents, 'Label', 'Set Label as');
  delete(get(setlabelas, 'children')); % remove existing list
  
  config = iView_private_config(instance, 'load');
  colors = iView_private_data_label(instance, 'colors');
  % add Default label item
  uimenu(setlabelas, 'Label', 'Default (no label)', 'Callback','iview(gcf, ''data_label_selection'',''Default'');');
  % then other labels
  for index=1:length(config.Labels)
  	label_index = mod(index, length(colors));
		label_index = colors{label_index};
  	uimenu(setlabelas, 'Label', config.Labels{index}, 'ForegroundColor', label_index, 'Callback',[ 'iview(gcf, ''data_label_selection'',''' config.Labels{index} ''');' ]);
  end

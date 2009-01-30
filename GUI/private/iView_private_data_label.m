function colors = iView_private_data_label(instance, action, label)
% iView_private_data_label label selected data sets in an iView instance
% action may be: add, remove, selection, colors, or Tags to iData/uicontrol
% if action=set, label indicates which label to set

% a list of not too violent colors to be used as labels
% from http://www.endprod.com/colors/ Web color definitions
colors{1} = [150  205  205]/256; % pale turquoise 3
colors{2} = [238  213  183]/256; % bisque 2
colors{3} = [152  245  255]/256; % cadet blue 2
colors{4} = [202  255  112]/256; % dark olive green 1
colors{5} = [193  255  193]/256; % dark sea green 2
colors{6} = [220  220  220]/256; % gainsboro
colors{7} = [240  230  140]/256; % khaki
colors{8} = [205  193  197]/256; % lavender blush 3
colors{9} = [205  201  165]/256; % lemon chiffon 3
colors{10} = [173  216  230]/256; % light blue
colors{11} = [240  128  128]/256; % light coral

config=iView_private_config(instance, 'load');

if ischar(action)

	switch(action)
	case 'add'
		  label = inputdlg(...
		    'Enter the name of the new Label to define:', ...
		    [ 'iView: Add a Label ' ]);
		  if ~isempty(label)
		    % add this Label if it is not defined in the list already
		    if isempty(strmatch(label{1}, config.Labels, 'exact'))
		      config.Labels{end+1} = label{1}; 
		    end
		    setappdata(0, 'iView_Config', config); % store new Labels
		    
		    % now update all Documents/Label menu in instances
		    instance_list=findall(0,'Tag','iView_instance');
		    for index=1:length(instance_list)
		      iView_private_documents(instance_list(index));
		    end
		  end
	case 'remove'
		  select = listdlg('PromptString', {'Select the Label(s)','to remove:'}, ...
		    'SelectionMode','multiple', ...
		    'OKstring', 'Remove', ...
		    'ListSize', [ 160 160 ], ...
		    'ListString', config.Labels, ...
		    'Name', 'iView: Remove a Label ');
		  instance_list=findall(0,'Tag','iView_instance');
		  labels = config.Labels;
		  for index=1:length(select)
		    label=labels{select(index)};
		    if strcmp(label,'Default'), continue; end % Default label can not be removed

		    % make sure the label is not used anymore
		    for instance_index = 1:length(instance_list)
		      Data = getappdata(instance_list(instance_index), 'Data');
		      Data_labels = get(Data,'Label');
		      if isempty(Data), continue; end
		      if ~iscellstr(Data_labels), Data_labels = cellstr(Data_labels); end
		      data_use_label = strmatch(label, Data_labels, 'exact');
		      if ~isempty(data_use_label)
		        % set label to none
		        iView_private_label(instance, data_use_label, '');
		      end
		    end
		    toremove = strmatch(label, config.Labels, 'exact');
		    config.Labels(toremove) = [];
		  end % for index
		  setappdata(0, 'iView_Config', config); % store new Labels
		  % now update all Documents/Label menu in instances
		  for index=1:length(instance_list)
		    iView_private_documents(instance_list(index));
		  end
	case 'selection'
		[selection, selectedIndex, selectedUI] = iView_private_selection(instance);
		iView_private_data_label(instance, selection, label);
	case 'colors'
  end
  return
else
  % use 'action' as a list of tags
  if isa(action, 'iData')
    tags = get(action,'Tag');
  elseif iscellstr(action)
  	tags = selection;
  end
  [selection, selectedIndex, selected] = iView_private_selection(instance, tags);
  if isempty(selection), return; end
  if strcmp(label,'Default'), label=''; end % Default label can not be removed
  selection = set(selection, 'Label', label);
  % update uicontrols
  for index=1:length(selection)
    this_data = ind2sub(selection, index);
    tooltip = iView_private_data_tooltip(this_data);
    icon = findobj(instance, 'Tag', this_data.Tag);
    set(icon, 'ToolTipString', tooltip);
  end % for
  % store these elements
  Data = getappdata(instance, 'Data');
  if length(Data) > 1
    Data(selectedIndex) = selection;
  else
    Data = selection;
  end
  setappdata(instance,'Data',Data);
  [hIcon,config,Data]=iView_private_icon(instance, 'check', Data);
end


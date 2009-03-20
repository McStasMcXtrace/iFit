function iView_private_data_rename(instance)
% iView_private_data_rename rename selected data sets in an iView instance
    [selection, selectedIndex, selectedUI] = iView_private_selection(instance);
  	if isempty(selection), return; end
  	this_name    = get(selection, 'Title');
  	this_tag     = get(selection, 'Tag');
  	if ischar(this_name), this_name = { this_name }; end
  	if ischar(this_tag),  this_tag  = { this_tag }; end
  	if length(selection) == 1, 
  		titl = 'a data set';
  	else
  		titl = [ num2str(length(selection)) ' data sets' ];
  	end
  	messg = strcat( 'Rename object:', this_tag );
  	options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='none';
  	newname = inputdlg(...
  		messg, ...
  		[ 'iView: Rename ' titl ' in window ' get(instance, 'Name') ], ...
  		1, this_name, options );
  	if ~isempty(newname)
  		% update icon text, tooltip, uicontectmenu and data sets
  		for index=1:length(selectedUI)
  			this_data = selection(index);
  			this_data.Title = newname{index};
  			if length(selection) > 1, selection(index) = this_data;
  			else selection = this_data; end
  			ud = get(selectedUI(index), 'UserData');
  			set(ud.uimenu_label, 'Label', newname{index});
  			set(selectedUI(index), 'String', newname{index});
  			tooltip = iView_private_data_tooltip(this_data);
				set(selectedUI(index), 'ToolTipString', tooltip);
  		end % for
  		% store these elements
  		Data = getappdata(instance, 'Data');
  		if length(Data) > 1
  			Data(selectedIndex) = selection;
  		else
  			Data = selection;
  		end
  		[hIcon,config,Data]=iView_private_icon(instance, 'check', Data);
  		setappdata(instance,'Data',Data);
    end

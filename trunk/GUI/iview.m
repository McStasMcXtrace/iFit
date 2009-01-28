function varargout=iview(varargin)
% IVIEW is a data browser, which allows to import any number of data sets, 
% create groups of data files, display/plot the data, apply mathematical 
% operations on single or multiple data sets, fit a model function to the data
% and export the data/plots.

% STORAGE of data
%   one iData array to hold data sets in the order of importation
%     each iData must have a reference to its corresponding checkbox/icon (handle)
%   uicontrol checkbox/icon, with reference to its corresponding iData (index in array and Tag)
%   one structure to hold configuration


% ARGUMENTS and ACTIONS
% new, load, save_data, resize open_data, close, exit, about, mouse_down, mouse_up, mouse_drag
% properties

varargout={};

% handling of arguments
if nargin==0, instance=[]; else instance=varargin{1}; end
if nargin<=1, action  =''; else action  =varargin{2}; end
if nargin<=2, object  =[]; else object  =varargin{3}; end
% check that instance is a single numeric figure handle
if ~isnumeric(instance)
  % shift arg to action
  if isempty(object), object=action; end
  action=instance;
  instance=[];
end
if ~ischar(action)
  % shift arg to object
  if isempty(object), object=action; end
  action='load';
end

% the interface is built in iView_private_create (menus, ...)
% create or raise interface
[instance, config]=iView_private_create(instance);

switch action
  case 'new'
    % open a new iView instance
    instance=figure;
    instance=iView_private_create(instance);
  case 'load'
    % open file selector to open data files
    set(instance, 'Pointer', 'watch');
    if isempty(object), object= iData(''); else object=iData(object); end
    [hIcon,config,object]=iView_private_icon(instance, 'load', object);
    iView_private_icon(instance, 'documents');
    set(instance, 'Pointer', 'arrow');
  case 'saveas_data'
    % open file selector to save data sets
    if ~isempty(object)
      set(instance, 'Pointer', 'watch');
      [hIcon,config,object]=iView_private_icon(instance, 'saveas', object);
      set(instance, 'Pointer', 'arrow');
    end
  case 'save_data'
    % open file selector to save data sets
    if ~isempty(object)
      set(instance, 'Pointer', 'watch');
      [hIcon,config,object]=iView_private_icon(instance, 'save', object);
      set(instance, 'Pointer', 'arrow');
    end
  case 'save_config'
    rmappdata(0,'iView_Config');
    iView_private_config(instance, 'save');
  case 'load_config'
  	config = iView_private_config(instance, 'load');
  	varargout{end+1} = config;
  case 'resize'
    % resize instance
    iView_private_icon(instance, 'resize',[]);
  case 'open_data'
    % open data set (plot it)
    [hIcon,config,object]=iView_private_icon(instance, 'open', object);
  case 'close_data'
    % close (delete) data sets
    [hIcon,config,object]=iView_private_icon(instance, 'delete', object);
  case 'close'
    % close this instance
    close(instance);
  case 'exit'
    if strcmp(config.ExitConfirm, 'yes') | config.ExitConfirm
      button = questdlg({'Do you really want to quit iView',' and close all windows ?'},'iView: Exit ?','Yes','No','Yes');
      if ~strcmp(button, 'Yes'), return; end
    end
    % close all instances
    % save configuration first
    iview(gcf, 'save_config');
    iLoad('','save config');
    % find all existing iView instances
    instance_list=findall(0,'Tag','iView_instance');
    close(instance_list);
    disp([ '% ' datestr(now) ' Exiting iView' ]);
  case 'about'
    % dialog box about iView
    helpstr = {[ 'iView ' config.Version ': This program enables to' ],...
	    '1- load data sets',...
	    '2- view/zoom/select data',...
	    '3- fit with any function or function set',...
	    '4- choose the fit method',...
	    '5- save your fit results and data',...
	    '6- configure many features at users choice',...
	    ' ',...
	    '*** iView comes with ABSOLUTELY NO WARRANTY',...
	    'This is free software, and you are welcome',...
	    'to redistribute it under certain conditions',...
	    'as specified in Licence files.'};
     helpdlg(helpstr,'iView: About');
  case 'contacts'
  	helpstr = {[ 'iView ' config.Version ],...
	    'Authors : ', ...
	    'E. Farhi <farhi@ill.fr>',...
	    'P.K.Willendrup <pkwi@risoe.dtu.dk>', ...
	    'Institut Laue-Langevin', ...
			'BP 156, 6, rue Jules Horowitz', ...
			'38042 Grenoble Cedex 9, France'};
	  helpdlg(helpstr,'iView: Contacts');
  case 'mouse_down'
    % handle mouse events (drag'n drop, move, select, ...)
    iView_private_mouse(instance, 'down', object);
  case 'mouse_drag'
    iView_private_mouse(instance, 'drag', object);
  case 'mouse_up'
    iView_private_mouse(instance, 'up', object);
  case 'select_all'
    iView_private_icon(instance, 'select_all', 1);
  case 'deselect_all'
    iView_private_icon(instance, 'deselect_all', 0);
  case 'properties'
    [hIcon,config,object]=iView_private_icon(instance, 'properties', object);
  case 'sort'
    items = {'Title (Name)','Size','Date','Label','Tag (unique ID)'};
    selection = listdlg('PromptString', {'Select sorting method to arrange data sets',[ 'in the iView window ' num2str(instance) ]}, ...
      'Name', 'iView: Sort data sets by...', ...
      'SelectionMode','single',...
      'OKString','Sort', 'ListSize', [160 150], ...
      'ListString', items);
    if isempty(selection), return; end
    Data = getappdata(gcf, 'Data');
    selection = items{selection};
    switch selection
    case 'Size'
      list = zeros(1,length(Data));
      for index=1:length(Data)
        list(index) = prod(size(Data(index)));
      end
    case 'Date'
      listd = get(Data,strtok(selection));
      for index=1:length(Data)
        list(index) = datenum(listd(index));
      end
    otherwise
      list = get(Data,strtok(selection));
    end
    [dummy, index] = sortrows(list(:));
    Data = ind2sub(Data, index);
    setappdata(instance, 'Data', Data);
    [hIcon,config,Data]=iView_private_icon(instance, 'check', Data);
    iView_private_icon(instance, 'documents', []);
  case 'rename_instance'
    this_name    = get(instance, 'Name');
		if this_name(1) == '*'
		  is_active = 1;
		  this_name(1) = '';
		else
			is_active = 0;
		end
  	newname = inputdlg(...
  		{ [ 'Enter the new name for the iView instance ' num2str(instance) ] }, ...
  		[ 'iView: Rename window ' num2str(instance) ], {this_name} );
  	if ~isempty(newname)
  		newname = newname{1};
  		if is_active, newname = [ '*' newname ]; end
  		set(instance, 'Name', newname);
  	end
  case 'rename_data'
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
  		% update uicontrols and data sets
  		for index=1:length(selectedUI)
  			this_data = ind2sub(selection, index);
  			set(selection(index), 'Title', newname{index});
  			ud = get(selectedUI(index), 'UserData');
  			set(ud.uimenu_label, 'Label', newname{index});
  			set(selectedUI(index), 'String', newname{index});
  			NL = sprintf('\n');
				tooltip = [ this_data.Title ' (' this_data.Label ')' NL this_data.Source NL 'Tag ' this_data.Tag '; Data size [' num2str(size(this_data)) ']' NL ];
				for ax_index = 0:ndims(this_data)
					[axisdef, lab] = getaxis(this_data, num2str(ax_index));
					lab=strtrim(lab);
					if ax_index==0, 
						if isempty(lab), tooltip = [ tooltip 'Signal' ]; else tooltip = [ tooltip lab ]; end
					elseif ~isempty(lab)
						tooltip = [ tooltip ' vs ' lab ]; 
					end
				end
				set(selectedUI(index), 'ToolTipString', tooltip);
  		end % for
  		% store these elements
  		Data = getappdata(gcf, 'Data');
  		if length(Data) > 1
  			Data(selectedIndex) = selection;
  		else
  			Data = selection;
  		end
  		setappdata(instance,'Data',Data);
  		[hIcon,config,Data]=iView_private_icon(instance, 'check', Data);
    end
  case 'label_data'
  	[selection, selectedIndex, selectedUI] = iView_private_selection(instance);
  	if isempty(selection), return; end
  	this_name    = get(selection, 'Title');
  	this_label   = get(selection, 'Label');
  	this_tag     = get(selection, 'Tag');
  	if ischar(this_name),  this_name = { this_name }; end
  	if ischar(this_label), this_label= { this_label }; end
  	if ischar(this_tag),   this_tag  = { this_tag }; end
  	if length(selection) == 1, 
  		titl = 'a data set';
  	else
  		titl = [ num2str(length(selection)) ' data sets' ];
  	end
  	messg = { [ 'Label ' titl '. WARNING: All selected objects' sprintf('\n') ...
  	'will be associated with the same label. Pre-defined labels as set' sprintf('\n') ...
  	'in File/Preferences will color data set icons.' ] };
  	options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='none';
  	newname = inputdlg(...
  		messg, ...
  		[ 'iView: Label ' titl ' in window ' get(instance, 'Name') ], ...
  		1, this_label(1), options );
  	if ~isempty(newname)
  		set(selection, 'Label', newname{1});
  		% update uicontrols
  		for index=1:length(selectedUI)
  			this_data = ind2sub(selection, index);
  			NL = sprintf('\n');
				tooltip = [ this_data.Title ' (' this_data.Label ')' NL this_data.Source NL 'Tag ' this_data.Tag '; Data size [' num2str(size(this_data)) ']' NL ];
				for ax_index = 0:ndims(this_data)
					[axisdef, lab] = getaxis(this_data, num2str(ax_index));
					lab=strtrim(lab);
					if ax_index==0, 
						if isempty(lab), tooltip = [ tooltip 'Signal' ]; else tooltip = [ tooltip lab ]; end
					elseif ~isempty(lab)
						tooltip = [ tooltip ' vs ' lab ]; 
					end
				end
				set(selectedUI(index), 'ToolTipString', tooltip);
  		end % for
  		% store these elements
  		Data = getappdata(gcf, 'Data');
  		if length(Data) > 1
  			Data(selectedIndex) = selection;
  		else
  			Data = selection;
  		end
  		setappdata(instance,'Data',Data);
  		[hIcon,config,Data]=iView_private_icon(instance, 'check', Data);
    end
  case 'cut'
    % copy
    % remove selected elements
  case 'copy'
    % copy selected data sets to internal clipboard buffer
    % set normal clipboard (for users to copy in other apps)
  case 'paste'
    % if external clipboard is different from the internal one, load files from clipboard('paste')
    % else load internal clipboard buffer into instance (if exists)
  case 'selection'
    varargout{end+1} = iView_private_selection(instance);
  case 'new_data'
  otherwise
    if ~isempty(action)
      disp([' Unknown action ' action ' in ' mfilename ]);
    end
end


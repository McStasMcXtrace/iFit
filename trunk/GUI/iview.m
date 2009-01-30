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
if ~strcmp(action, 'instance_close')
	[instance, config]=iView_private_create(instance);
end

switch action

% ACTIONS on instance
  case 'instance_new'
    % open a new iView instance
    instance=figure;
    instance=iView_private_create(instance);  
  case 'instance_resize'
    % resize instance and re-arrange icons
    iView_private_icon(instance, 'resize',[]);
  case 'instance_sort'
    % sort data sets in instance
    iView_private_sort(instance, object);
  case 'instance_rename'
    this_name    = get(instance, 'Name');
		if this_name(1) == '*'
		  is_active = 1;
		  this_name(1) = '';
		else
			is_active = 0;
		end
  	newname = inputdlg(...
  		{ [ 'Enter the new name for the iView instance ' num2str(instance) ] }, ...
  		[ 'iView: Rename window ' num2str(instance) ], 1, {this_name} );
  	if ~isempty(newname)
  		newname = newname{1};
  		if is_active, newname = [ '*' newname ]; end
  		set(instance, 'Name', newname);
  	end
  case 'instance_close'
    % find all existing iView instances
    instance_list=findall(0,'Tag','iView_instance');
    if length(instance_list) == 1 % last instance to be closed
      disp([ '% ' datestr(now) ' Exiting iView' ]);
    	% save configuration first
    	iview(gcf, 'config_save');
    	rmappdata(0,'iView_Config');
    end
    % close this instance
    delete(instance);
  case 'exit'
    if strcmp(config.ExitConfirm, 'yes') | config.ExitConfirm
      button = questdlg({'Do you really want to quit iView',' and close all windows ?'},'iView: Exit ?','Yes','No','Yes');
      if ~strcmp(button, 'Yes'), return; end
    end
    % close all instances
    % save configuration first
    iview(gcf, 'config_save');
    % find all existing iView instances
    instance_list=findall(0,'Tag','iView_instance');
    close(instance_list);
% ACTIONS on data    
  case 'data_load'
    % open file selector to open data files
    set(instance, 'Pointer', 'watch');
    if isempty(object), object= iData(''); else object=iData(object); end
    % loading is called from iView_private_icon with iData(char)/uigetfiles
    [hIcon,config,object]=iView_private_icon(instance, 'load', object);
    iView_private_documents(instance);
    set(instance, 'Pointer', 'arrow');
  case 'data_saveas'
    % open file selector to save data sets
    if ~isempty(object)
      set(instance, 'Pointer', 'watch');
      [hIcon,config,object]=iView_private_icon(instance, 'saveas', object);
      set(instance, 'Pointer', 'arrow');
    end
  case 'data_save'
    % open file selector to save data sets
    if ~isempty(object)
      set(instance, 'Pointer', 'watch');
      [hIcon,config,object]=iView_private_icon(instance, 'save', object);
      set(instance, 'Pointer', 'arrow');
    end  
  case 'data_open'
    % open data set (plot it)
    [hIcon,config,object]=iView_private_icon(instance, 'open', object);
  case 'data_close'
    % close (delete) data sets
    [hIcon,config,object]=iView_private_icon(instance, 'delete', object);
  case 'data_select_all'
    iView_private_icon(instance, 'select_all', 1);
  case 'data_deselect_all'
    iView_private_icon(instance, 'deselect_all', 0);
  case 'data_properties'
    [hIcon,config,object]=iView_private_icon(instance, 'properties', object);
  case 'data_rename'
  	iView_private_data_rename(instance);
  case 'data_label' % add/remove
  	iView_private_data_label(instance, object);
  case 'data_label_selection' % add/remove
  	iView_private_data_label(instance, 'selection', object);
  case 'data_cut'
    % copy
    % remove selected elements
  case 'data_copy'
    % copy selected data sets to internal clipboard buffer
    % set normal clipboard (for users to copy in other apps)
  case 'data_paste'
    % if external clipboard is different from the internal one, load files from clipboard('paste')
    % else load internal clipboard buffer into instance (if exists)
  case 'selection'
    varargout{end+1} = iView_private_selection(instance);
  case 'data_new'
  
% ACTIONS on config    
  case 'config_save'
    if isappdata(0,'iView_Config')
      rmappdata(0,'iView_Config');
    end
    iView_private_config(instance, 'save');
    iLoad_config=iLoad('', 'load config');
    iLoad_config.UseSystemDialogs = config.UseSystemDialogs;
    iLoad(iLoad_config, 'save config');
  case 'config_load'
  	config = iView_private_config(instance, 'load');
  	varargout{end+1} = config;
 
% ACTIONS with mouse	  
  case 'mouse_down'
    % handle mouse events (drag'n drop, move, select, ...)
    iView_private_mouse(instance, 'down', object);
  case 'mouse_drag'
    iView_private_mouse(instance, 'drag', object);
  case 'mouse_up'
    iView_private_mouse(instance, 'up', object);
    
% OTHER actions
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

  otherwise
    if ~isempty(action)
      disp([' Unknown action ' action ' in ' mfilename ]);
    end
end


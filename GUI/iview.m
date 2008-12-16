function iview(varargin)
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
instance=iView_private_create(instance);

switch action
  case 'new'
    % open a new iView instance
    instance=figure;
    instance=iView_private_create(instance);
  case 'load'
    % open file selector to open data files
    set(instance, 'Pointer', 'watch');
    if isempty(object), object= iData(''); else object=iData(object); end
    iView_private_icon(instance, 'load', object);
    iView_private_icon(instance, 'documents', []);
    set(instance, 'Pointer', 'arrow');
  case 'save_data'
    % open file selector to save data sets
    if ~isempty(object)
      set(instance, 'Pointer', 'watch');
      iView_private_icon(instance, 'saveas', object);
      set(instance, 'Pointer', 'arrow');
    end
  case 'save_config'
    iView_private_config(instance, 'save');
  case 'resize'
    % resize instance
    iView_private_icon(instance, 'resize', []);
  case 'open_data'
    % open data set (plot it)
    iView_private_icon(instance, 'open', object);
  case 'close_data'
    % close (delete) data se
    iView_private_icon(instance, 'delete', object);
  case 'close'
    % close this instance
    close(instance);
  case 'exit'
    % close all instances
    % save configuration first
    iview(gcf, 'save_config');
    % find all existing iView instances
    iLoad('','save config');
    instance_list=findall(0,'Tag','iView_instance');
    close(instance_list);
    disp([ '% ' datestr(now) ' Exiting iView' ]);
  case 'about'
    % dialog box about iView
    helpstr = {'iView 0.1: This program enables to',...
	  '1- load data sets',...
	  '2- view/zoom/select data',...
	  '3- fit with any function or function set',...
	  '4- choose the fit method',...
	  '5- save your fit results and data',...
	  '6- configure many features at users choice',...
	  ' ',...
	  'Authors : E. Farhi <farhi@ill.fr>',...
	  '*** iView comes with ABSOLUTELY NO WARRANTY',...
	  'This is free software, and you are welcome',...
	  'to redistribute it under certain conditions',...
	  'as specified in Licence files.'};
     helpdlg(helpstr,'iView About');
  case 'mouse_down'
    % handle mouse events (drag'n drop, move, select, ...)
    iView_private_mouse(instance, 'down', object);
  case 'mouse_drag'
    iView_private_mouse(instance, 'drag', object);
  case 'mouse_up'
    iView_private_mouse(instance, 'up', object);
  otherwise
    if ~isempty(action)
      disp([' Unknown action ' action ' in ' mfilename ]);
    end
end


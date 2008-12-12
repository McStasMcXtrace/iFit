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

% handling of arguments
if nargin==0, instance=[]; else instance=varargin{1}; end
if nargin<=1, action  =''; else action  =varargin{2}; end
if nargin<=2, object  =[]; else object  =varargin{3}; end
% check that instance is a single numeric figure handle
if ~isnumeric(instance)
  % shift arg to action
  if isempty(action), action=instance; end
  instance=[];
end
if ~ischar(action)
  % shift arg to object
  if isempty(object), object=action; end
  action='load';
end

% create or raise interface
instance=iView_private_create(instance);

switch action
case 'new'
  instance=figure;
  instance=iView_private_create(instance);
  return
case 'load'
  if isempty(object), object= iData(''); end
  iView_private_icon(instance, 'load', object)
case 'resize'
  iView_private_icon(instance, 'resize', []);
  return
otherwise
  if ~isempty(action)
    disp([' Unknown action ' action ' in ' mfilename ]);
  end
end
% iview actions:
%   build         build empty interface
%   load_config   read INI file
%   save_config   write INI file (config only, not the data)
%   load          send argument to iData, import and create/update icons
%   save          save config and all data sets into INI/M script
%   delete        remove data set and icon
%   print         print interface or selected elements
%   select        toggle icon selection state
%   update        check all data sets and update icons/menus
%   align         re-arrange all icons and align them on a grid
%   exit          save configuration (if needed), clear memory and destroy interface
%   copy          identify data set and store it into clipboard
%   paste         paste data set
%   duplicate     copy+paste
%   edit          edit iView configuration or data sets (if selection is active)

% INTERFACE



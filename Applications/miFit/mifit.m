function varargout = mifit(varargin)


% notes:
% saveas: display dialogue to select what to export:
% * data sets
% * parameters
% * fit model

% plot: plot model/dataset/parameters ? -> preferences
% will add toolbar instead of uicontrols for fast processing

    varargout = {};
    % look if the main window is opened
    fig = mifit_fig();
    if isempty(fig) || ~ishandle(fig)
        fig = feval([ mfilename '_OpeningFcn' ]);
    end
    out = fig;

    if ~isempty(varargin)
        if ischar(varargin{1}) && isempty(dir(varargin{1})) % a function/action to call ?
          % callback with varargin{1} == 'action'
          action = varargin{1};
          if strcmpi(action,'identify'), return; end
          try
            feval([ 'mifit_' action ], varargin{2:end});
          catch ME
            mifit_disp([ mfilename ': Unknown action "' action '"' ])
            rethrow(ME)
          end
        elseif ishandle(varargin{1})    % a CallBack function
          h = varargin{1}; d = [];
          if strcmpi(get(h,'Type'),'uitable')
              % get data from a UITable
              try
                  d = edit(iData, h);
              catch
                  d = iData(get(h, 'Data'));
                  d.Title = get(h, 'Tag');
              end
              % now push 'd' into the Stack, except for 'empty' Table
              if numel(d) == 1 && numel(find(d == 0)) == 25 % zeros(5)
                d = [];
              end
          else
            d = iData(h);
          end
          if ~isempty(d)
            mifit_disp([ mfilename ': Import into List:' ])
            mifit_disp(char(d))
            mifit_List_Data_push(d);
          end
        else
          d = iData(varargin{:});
          % now push 'd' into the Stack
          if ~isempty(d)
            mifit_disp([ mfilename ': Import into List:' ])
            mifit_disp(char(d))
            mifit_List_Data_push(d);
          end
        end
    end
    if nargout >= 1
        varargout{end+1} = out;
    end
    
% -------------------------------------------------------------------------
% Main window ID and creation
% -------------------------------------------------------------------------

function f=mifit_fig(tag)
% search for a given Tag in Application or main Figure if ommitted.
  persistent fig handles
  
  if ~ishandle(fig), fig=[]; end
  if isempty(fig)
    fig = findall(0, 'Tag','miFit');
    if length(fig) > 1, delete(fig(2:end)); end % unique instance
    handles = [];
  end

  if nargin == 0
    f=fig;
  else
    if ~isfield(handles, tag) || ~ishandle(handles.(tag))
      handles.(tag) = findobj(fig, 'Tag', tag);
    end
    f = handles.(tag);
  end

% --- Executes just before mifit is made visible.
function fig = mifit_OpeningFcn
% This function creates the main window and returns its ID

fig = mifit_fig();
if isempty(fig) || ~ishandle(fig)
    % create the main figure
    mifit_disp([ mfilename ': Welcome to miFit !' ])
    mifit_disp(datestr(now))
    fig = openfig(mfilename);
    
    % get Models/Optimizers menu handles
    hmodels = mifit_fig('Menu_Models');
    hoptim  = mifit_fig('Menu_Optimizers');
    
    % Display welcome dialog during menu build
    h = mifit_Tools_About(fig);
    
    % get the list of Models and Optimizers
    [optimizers,functions] = fits(iFunc);
    
    % fill Models menu
    if any(~isempty(functions)) && all(isa(functions, 'iFunc'))
        mifit_disp([ mfilename ': Initializing Models... User Models should be in the local directory.' ]);
        separator = 1;
        for f=functions
            % each Model is an iFunc object. These should be stored in the
            % Models items 'UserData'
            if ~isempty(f) && ~isempty(f.Name)
                if separator
                  uimenu(hmodels, 'Label', f.Name, 'UserData', f, 'Separator','on');
                  separator = 0;
                else
                  uimenu(hmodels, 'Label', f.Name, 'UserData', f);
                end
            end
        end
    end
    
    % fill Optimizers menu
    if ~isempty(optimizers) && iscell(optimizers)
        mifit_disp([ mfilename ': Initializing Optimizers...' ]);
        for f=optimizers
            % each optimizer is given with its function name. We request
            % 'defaults' and display its name
            
            o=feval(f{1},'defaults');
            if isfield(o, 'algorithm') && ~isempty(o.algorithm)
                algorithm = o.algorithm;
            else
                algorithm = f{1};
            end
            if ~isempty(algorithm)
                uimenu(hoptim, 'Label', algorithm, 'UserData', f{1});
            end
        end
    end
    
    % create the AppData Data Stack
    Data = [];
    setappdata(fig, 'Data', Data);
    setappdata(fig, 'History', {});
    
    % Load the save configuration: Data sets and Model Parameters
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file))
      try
        d = load(file);
        mifit_disp([ mfilename ': Loading previous session Data sets from ' file ]);
        if isfield(d, 'Data')
          Data = d.Data;
          mifit_List_Data_push(Data);
        end
      end
    end

    % close welcome image
    delete(h);
end

% -------------------------------------------------------------------------
% Callbacks
% -------------------------------------------------------------------------

% File menu ********************************************************************
function mifit_File_New(handle)
% File/New menu item
  d = iData(zeros(5)); % create an empty Data set;
  mifit_disp([ mfilename ': Editing an empty Data set. Close the window to retrieve its content' ]);
  mifit_disp( '    as a new Data set into miFit. Use the Contextual menu.');
  handle = edit(d, 'editable');
  % set a DeletedFcn so that the content can be retrieved into miFit when
  % closing.
  set(handle, 'DeleteFcn', @mifit);
  
function mifit_File_Open(handle)
  d = iData('');  % open file selector, and import files
  % push that data onto the List
  mifit_List_Data_push(d);
  
function mifit_File_Save(hObject)
% save Data sets and Model parameters into a mifit.mat file
  fig = mifit_fig;
  Data = getappdata(fig, 'Data');
  if isempty(Data), return; end
  
  file = fullfile(prefdir, [ mfilename '.mat' ]);
  mifit_disp([ mfilename ': Saving Data sets into ' file ]);
  builtin('save', file, 'Data');
  
function mifit_File_Saveas(hObject)
% save the application configuration into specified file
% Data sets and Model parameters into a mifit.mat file
  fig = mifit_fig;
  Data = getappdata(fig, 'Data');
  if isempty(Data), return; end
 
  filterspec = { '*.mat','MAT-files (*.mat)'; ...
                 '*.html','Web HTML document' };
  [filename, pathname, filterindex] = uiputfile('Save miFit workspaces as', [ mfilename '.mat' ]);
  if isequal(filename,0) || isequal(pathname,0)
    return
  end
  file = fullfile(pathname, filename);
  mifit_disp([ mfilename ': Saving Data sets into ' file ]);
  switch filterindex
  case 1
    builtin('save', file, 'Data');
  case 2
    % write an HTML report
    mifit_File_Saveas_HTML(file); % display list of data sets, plots, models, model parameters
  end
  
function mifit_File_Print(hObject)
% print the interface. Not very nice. can we think of something better ?
% perhaps we can generate an HTML report in Saveas HTML ?
  fig = mifit_fig;
  printdlg(fig);
  % alternative: File_Saveas_HTML in tmpfile, then open that file for printing.
  disp('TODO: save all Data sets as HTML with model and parameters')
  disp('then open it with web for printing');
  
function mifit_File_Preferences(hObject)
% open Preferences dialogue
% set directories to search for Models
% set FontSize (and update all Fonts in figure)
% set Save on exit
% save Preferences on dialogue close
  fig = mifit_fig;
  promt = {'Font size','Save Data sets on Exit'};
  defaultanswer = { get(fig, 'FontSize'), 'yes' };
  name  = [ mfilename ': Preferences' ];
  options.Resize='on';
  options.WindowStyle='normal';
  disp('TODO: mifit_File_Preferences')

function mifit_File_Exit(hObject)
% Quit and Save Data
  mifit_File_Save;
  mifit_disp([ mfilename ': Exiting miFit. Bye bye.' ])
  mifit_disp(datestr(now))
  delete(mifit_fig);
  
% Edit menu ********************************************************************

function mifit_Edit_Undo(hObject)
% set the Data stack to the previous state from History
  mifit_History_pull();
  mifit_List_Data_UpdateStrings();

function mifit_Edit_Cut(hObject)
% get the selected indices in the List, copy these elements to the clipboard
% and delete the elements. Update the History (in Delete).
  mifit_Copy(hObject);
  mifit_Delete(hObject);

function mifit_Edit_Copy(hObject)
% get the selected indices in the List, copy these elements to the clipboard
  

function mifit_Edit_Paste(hObject)
% append/copy the data sets from the clipboard to the end of the list
% clipboard can be a file name, an iData Tag/ID or ID in the Stack
% for iData ID, make a copy of the objects
% Update the History
  disp('TODO: mifit_Edit_Paste: ...')
  d=[];
  mifit_List_Data_push(d);
  mifit_History_push;

function mifit_Edit_Duplicate(hObject)
% copy and paste selected. Update the history (in Paste).
  mifit_Copy(hObject);
  mifit_Paste(hObject);
  
function mifit_Edit_Select_All(hObject)
% set the List selected values to all ones
  hObject = mifit_fig('List_Data_Files');
  items   = get(hObject,'String');
  index_selected = get(hObject,'Value');
  if numel(index_selected) == numel(items), index_selected = [];
  else index_selected=1:numel(items); end
  set(hObject,'Value', index_selected);

function mifit_Edit_Delete(hObject)
% delete selected
% Update the History
  fig = mifit_fig;
  hObject        = mifit_fig('List_Data_Files');
  index_selected = get(hObject,'Value');
  Data = getappdata(fig, 'Data');
  if numel(Data) > 1
    Data(index_selected) = [];
  else
    Data = [];
  end
  list           = get(hObject,'String');
  list(index_selected) = [];
  set(hObject,'Value',[]);
  set(hObject,'String',list);
  setappdata(fig, 'Data', Data);
  
  mifit_History_push;

function mifit_History_pull
% get the last History element and deletes it
  fig = mifit_fig;
  History = getappdata(fig, 'History');
  if isempty(History), return; end
  Data         = History{end};
  History(end) = [];
  setappdata(fig, 'History', History);
  setappdata(fig, 'Data',    Data);
  

function mifit_History_push
% append current Data into the History
  fig = mifit_fig;
  History = getappdata(fig, 'History');
  Data    = getappdata(fig, 'Data');
  History{end+1} = Data;
  if numel(History) > 10, History(1:(end-9)) = []; end
  setappdata(fig, 'History', History);
  
% Data menu ********************************************************************

% Data_Properties ?
% re-assign data set signal, axes, ... to aliases/new ones
% display statistics

% Models and Optimizers menu ***************************************************

% set optimizer configuration -> contextual dialogue
  
% Tools menu *******************************************************************

function h=mifit_Tools_About(fig)
% display the About dialogue. The handle ID is in adddata(gcf, 'handle_About')
  if nargin ==0, fig=''; end
  icon = fullfile(ifitpath,'Docs','images','ILL-web-jpeg.jpg');
  
  % Display About dialog
  t = [ sprintf('Welcome to miFit, a GUI to iFit.\n ') version(iData,2) sprintf('.\n Visit <http://ifit.mccode.org>') ];
  if isempty(dir(icon))
    h = msgbox(t,'miFit: About','help');
  else
    h = msgbox(t,'miFit: About','custom', imread(icon));
  end
  if ~isempty(fig)
    setappdata(fig, 'handle_About', h);
  end
  
% List Data and Stack management ***********************************************

function mifit_List_Data_push(d)
% put a new data set at the end of the stack
  if isempty(d), return; end
  fig = mifit_fig;

  % update AppData Stack
  if numel(d) > 1, d = d(:); end
  Data = getappdata(fig, 'Data');
  Data = [ Data ; d ];  % a column of iData set
  setappdata(fig, 'Data', Data);
  
  % update the List labels by appending the Name at the end
  hObject        = mifit_fig('List_Data_Files');
  list           = get(hObject,'String');
  list0          = numel(list);
  index_selected = get(hObject,'Value');
  if max(index_selected) > numel(list), index_selected = []; end
  for index=1:numel(d)
      index_selected(end+1) = list0+index;
      list{end+1} = char(d(index));
  end
  set(hObject,'String', list, 'Value', index_selected);
  
  % Update the History with the new stack and Date
  mifit_History_push;
  
function d=mifit_List_Data_pull(hObject)
% get the selected Data List
% return the selected objects
  if nargin == 0 || isempty(hObject)
      hObject = mifit_fig('List_Data_Files');
  end
  d = [];

  index_selected = get(hObject,'Value');
  if isempty(index_selected), return; end
  
  fig = mifit_fig;
  d   = getappdata(fig, 'Data');
  if numel(d) > 1
      d = d(index_selected);
  end
  
function mifit_List_Data_UpdateStrings
  % update the List labels
  fig = mifit_fig;
  
  Data = getappdata(fig, 'Data');
  
  hObject        = mifit_fig('List_Data_Files');
  list           = {};
  for index=1:numel(Data)
      list{end+1} = char(Data(index));
  end
  set(hObject,'String', list, 'Value', []);

function mifit_disp(message)
  % display message, and can log it into a File or Log window
  disp(message)
  disp('TODO: mifit_disp: log into file, set in Preferences')

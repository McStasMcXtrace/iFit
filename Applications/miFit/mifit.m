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
    fig = mifit_fig;
    if isempty(fig) || ~ishandle(fig)
        fig = feval([ mfilename '_OpeningFcn' ]);
    end
    out = fig;

    if nargin == 1 && isa(varargin{1}, 'iFunc')
        % import model: add it to the Models menu, if not already there
    elseif nargin == 1 && isa(varargin{1}, 'iData')
        % import iData object(s) into the List
        disp([ mfilename ': Import into List:' ])
        disp(char(varargin{1}))
        mifit_List_Data_push(varargin{1});
    end

    if ~isempty(varargin)
        if ischar(varargin{1}) && isempty(dir(varargin{1})) % a function/action to call ?
            % callback with varargin{1} == 'action'
            action = varargin{1};
            if strcmpi(action,'identify'), return; end
            try
              feval([ 'mifit_' action ], varargin{2:end});
            catch ME
              disp([ mfilename ': Unknown action "' action '"' ])
              rethrow(ME)
            end
        elseif ishandle(varargin{1})    % a CallBack function
            h = varargin{1};
            if strcmpi(get(h,'Type'),'uitable')
                % get data from a UITable
                try
                    d =edit(iData, h);
                catch
                    d = iData(get(h, 'Data'));
                    d.Title = get(h, 'Tag');
                end
                % now push 'd' into the Stack
                mifit_List_Data_push(d);
            end
        else
          d = iData(varargin{:})
          % now push 'd' into the Stack
          mifit_List_Data_push(d);
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

fig = mifit_fig;
if isempty(fig) || ~ishandle(fig)
    % create the main figure
    disp([ mfilename ': Welcome to miFit !' ])
    disp(datestr(now))
    fig = openfig(mfilename);
    
    % test if any of the Models/Optimizers menu is empty
    % if not, return (no need to build them)
    hmodels = mifit_fig('Menu_Models');
    hoptim  = mifit_fig('Menu_Optimizers');
    if isempty(hmodels) || isempty(hoptim) ...
            || ~isempty([ get(hmodels,'Children') ; get(hoptim,'Children') ]), 
      return; 
    end
    
    % Display welcome dialog during menu build
    h = mifit_Tools_About(fig);
    
    % get the list of Models and Optimizers
    [optimizers,functions] = fits(iFunc);
    
    % fill Models menu
    if any(~isempty(functions)) && all(isa(functions, 'iFunc'))
        disp([ mfilename ': Initializing Models...' ]);
        separator = 1;
        % first add the 'Add new Model': from File or ifitmakefunc
        uimenu(hmodels, 'Label','Add new Model...');
        uimenu(hmodels, 'Label','Plot Model...');
        uimenu(hmodels, 'Label','View Model Parameters...', 'Separator','on');
        uimenu(hmodels, 'Label','Plot Model Parameters...');
        uimenu(hmodels, 'Label','Export Model Parameters...');
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
        disp([ mfilename ': Initializing Optimizers...' ]);
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
    
    % Load the save configuration: Data sets and Model Parameters
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file))
      try
        d = load(file)
        disp([ mfilename ': Loading previous session Data sets from ' file ]);
        Data = d.Data;
        clear d
      end
    end
    if ~isa(Data, 'iData'), Data=[]; end
    setappdata(fig, 'Data', Data);
    
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
  disp([ mfilename ': Editing an empty Data set. Close the window to retrieve its content' ]);
  disp( '    as a new Data set into miFit. Use the Contextual menu.');
  handle = edit(d, 'editable');
  % set a DeletedFcn so that the content can be retrieved into miFit when
  % closing.
  set(handle, 'DeleteFcn', @mifit);
  
function mifit_File_Open(handle)
  d = iData('');
  % push that data onto the List
  mifit_List_Data_push(d);
  
function mifit_File_Save(hObject)
% save Data sets and Model parameters into a mifit.mat file
  fig = mifit_fig;
  Data = getappdata(fig, 'Data');
  if isempty(Data), return; end
  
  file = fullfile(prefdir, [ mfilename '.mat' ]);
  disp([ mfilename ': Saving Data sets into ' file ]);
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
  disp([ mfilename ': Saving Data sets into ' file ]);
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
  
function mifit_File_Preferences(hObject)
% open Preferences dialogue
% set direcories to search for Models
% set FontSize (and update all Fonts in figure)
% save Preferences on dialogue close

function mifit_File_Exit(hObject)
% Quit and Save Data
  mifit_File_Save;
  disp([ mfilename ': Exiting miFit. Bye bye.' ])
  disp(datestr(now))
  
% Edit menu ********************************************************************

function mifit_Edit_Undo(hObject)
% set the Data stack to the previous state from History

function mifit_Edit_Cut(hObject)
% get the selected indices in the List, copy these elements to the clipboard
% and delete the elements. Update the History (in Delete).
  mifit_Copy(hObject);
  mifit_Delete(hObject);

function mifit_Edit_Copy(hObject)
% get the selected indices in the List, copy these elements to the clipboard

function mifit_Edit_Paste(hObject)
% append/copy the data sets from the clipboard to the end of the list
% clipboard can be a file name, an iData Tag/ID or the ID in the Stack
% Update the History

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

function mifit_Delete(hObject)
% delete selected
% Update the History

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
  fig = mifit_fig;

  % update AppData Stack
  if numel(d) > 1, d = d(:); end
  Data = getappdata(fig, 'Data');
  Data = [ Data ; d ];  % a columns of iData set
  setappdata(fig, 'Data', Data);
  
  % update the List labels by appending the Name at the end
  hObject        = mifit_fig('List_Data_Files');
  list           = get(hObject,'String');
  index_selected = get(hObject,'Value');
  for index=1:numel(d)
      list{end+1} = char(d(index));
      index_selected(end+1) = numel(list)+index;
  end
  set(hObject,'String', list, 'Value', index_selected);
  
  % Update the History with the new stack and Date
  
  
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
  

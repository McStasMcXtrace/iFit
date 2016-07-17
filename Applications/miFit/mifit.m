function varargout = mifit(varargin)


% notes:
% saveas: display dialogue to select what to export:
% * data sets
% * parameters
% * fit model

% plot: plot model/dataset/parameters ? -> preferences
% will add toolbar instead of uicontrols for fast processing
    persistent config
    
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
          if any(strcmpi(action,{'pull','data'}))
            out = getappdata(fig, 'Data');
          else
            try
              feval([ 'mifit_' action ], varargin{2:end});
            catch ME
              mifit_disp([ 'Unknown action "' action '"' ])
              rethrow(ME)
            end
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
            mifit_List_Data_push(d);
          end
        else
          d = iData(varargin{:});
          % now push 'd' into the Stack
          if ~isempty(d)
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

% --- Executes just before mifit is made visible -------------------------------
function fig = mifit_OpeningFcn
% This function creates the main window and returns its ID

fig = mifit_fig();
if isempty(fig) || ~ishandle(fig)
    % create the main figure
    mifit_disp('Welcome to miFit ! ********************************************')
    fig = openfig(mfilename);
    
    % get Models/Optimizers menu handles
    hmodels = mifit_fig('Menu_Model');
    hoptim  = mifit_fig('Menu_Optimizers');
    
    % Display welcome dialog during menu build
    h = mifit_Tools_About(fig);
    
    % get the list of Models and Optimizers
    [optimizers,functions] = fits(iFunc);
    
    % fill Models menu
    if any(~isempty(functions)) && all(isa(functions, 'iFunc'))
        mifit_disp([ 'Initializing Models... User Models should be in the local directory.' ]);
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
        mifit_disp([ 'Initializing Optimizers...' ]);
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
    setappdata(fig, 'Data',    []);
    setappdata(fig, 'History', {});
    
    % load Preferences
    mifit_Load_Preferences;
    mifit_Apply_Preferences;
    
    % Load the previous Data sets and Model Parameters
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file))
      try
        d = load(file);
        if isfield(d, 'Data')
          mifit_disp([ 'Loading Data sets from ' file ]);
          mifit_List_Data_push(d.Data);
        end
      end
    end
    file = fullfile(prefdir, [ mfilename '.log' ]);
    mifit_disp([ 'Log file is ' file ]);

    % close welcome image
    delete(h);
end

% Preferences I/O --------------------------------------------------------------
function config = mifit_Load_Preferences
  file = fullfile(prefdir, [ mfilename '.ini' ]);
  content = ''; config = '';
  if ~isempty(dir(file))
    try
      content = fileread(file);
      evalc(content);% this should make a 'config' variable
      mifit_disp([ 'Loading Preferences from ' file ]);
    end
  end
  if isempty(config)
    % default configuration
    config.FontSize         = 12;
    config.Save_Data_On_Exit= 'yes';
  end
  setappdata(mifit_fig, 'Preferences', config);

function mifit_Save_Preferences(config)
  filename = fullfile(prefdir, [ mfilename '.ini' ]);
  NL = sprintf('\n');
  description = 'miFit interface to iFit';
    str = [ '% miFit configuration script file: ' description NL ...
          '%' NL ...
          '% Matlab ' version ' m-file ' filename NL ...
          '% generated automatically on ' datestr(now) ' with ifit.mccode.org ' mfilename NL...
          class2str('config', config) ];
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    warning([ datestr(now) ': Error opening file ' filename ' to save ' description ' configuration.' ]);
    filename = [];
  else
    fprintf(fid, '%s', str);
    fclose(fid);
    mifit_disp([ 'Saving Preferences into ' filename ]);
  end

% -------------------------------------------------------------------------
% Callbacks
% -------------------------------------------------------------------------
% usual callback arguments are : Callback(gcbo,eventdata)

% File menu ********************************************************************
function mifit_File_New(handle)
% File/New menu item
  d = iData(zeros(5)); % create an empty Data set;
  mifit_disp([ 'Editing an empty Data set. ' ...
    'Close the window to retrieve its content as a new Data set into miFit. ' ...
    'Use the Contextual menu for Copy/Paste/Resize.' ]);
  handle = edit(d, 'editable');
  % set a DeletedFcn so that the content can be retrieved into miFit when
  % closing.
  set(handle, 'DeleteFcn', @mifit);
  
function mifit_File_Open(handle)
  d = iData('');  % open file selector, and import files
  % push that data onto the List
  mifit_List_Data_push(d);
  
function mifit_File_Save(varargin)
% save Data sets and Model parameters into a mifit.mat file
  fig = mifit_fig;
  Data = getappdata(fig, 'Data');
  if isempty(Data), return; end
  
  file = fullfile(prefdir, [ mfilename '.mat' ]);
  mifit_disp([ 'Saving Data sets into ' file ]);
  builtin('save', file, 'Data');
  
function mifit_File_Saveas(varargin)
% save the application configuration into specified file
% Data sets and Model parameters into a mifit.mat file
  fig  = mifit_fig;
  Data = getappdata(fig, 'Data');
  if isempty(Data), return; end
 
  filterspec = { '*.mat','MAT-files (*.mat)'};
  [filename, pathname] = uiputfile(filterspec, 'Save All miFit Data sets as', [ mfilename '.mat' ]);
  if isequal(filename,0) || isequal(pathname,0)
    return
  end
  file = fullfile(pathname, filename);
  mifit_disp([ 'Saving Data sets into ' file ]);
  builtin('save', file, 'Data');
  
function mifit_File_Print(varargin)
% print the interface. Not very nice. can we think of something better ?
% perhaps we can generate an HTML report in Saveas HTML ?
  fig = mifit_fig;
  printdlg(fig);
  % alternative: File_Saveas_HTML in tmpfile, then open that file for printing.
  disp('TODO: save all Data sets as HTML with model and parameters')
  disp('then open it with web for printing');
  
function mifit_File_Preferences(varargin)
% open Preferences dialogue
% set directories to search for Models
% set FontSize (and update all Fonts in figure)
% set Save on exit
% save Preferences on dialogue close
  fig = mifit_fig;
  config = getappdata(mifit_fig, 'Preferences');
  prompt = {'Font size [10-36]','Save Data sets on Exit [yes/no]'};
  defaultanswer = { num2str(config.FontSize), config.Save_Data_On_Exit };
  name  = [ mfilename ': Preferences' ];
  options.Resize='on';
  options.WindowStyle='normal';
  answer=inputdlg(prompt,name,1,defaultanswer,options);
  if isempty(answer), return; end
  
  % set new Preferences
  answer{1} = str2double(answer{1});
  if isfinite(answer{1}) 
    config.FontSize=min(max(answer{1}, 10),36); end
  if any(strcmp(answer{2}, {'yes','no'})) 
    config.Save_Data_On_Exit = answer{2}; end
  setappdata(mifit_fig, 'Preferences', config);
  mifit_Apply_Preferences;
  mifit_Save_Preferences(config);
  
function mifit_Apply_Preferences
  fig = mifit_fig;
  config = getappdata(fig, 'Preferences');
  % change all Font Size
  h=[ findobj(fig, 'Type','uicontrol') ; ...
      findobj(fig, 'Type','axes') ; findobj(fig, 'Type','text') ; ...
      findobj(fig, 'Type','uipanel') findobj(fig, 'Type','uitable') ];
  set(h, 'FontSize', config.FontSize);
  % for uimenu, we use tip given by Y Altman 
  % <https://fr.mathworks.com/matlabcentral/newsreader/view_thread/148095>
  h = findobj(fig, 'Type','uimenu');

function mifit_File_Exit(varargin)
% Quit and Save Data
  config = getappdata(mifit_fig, 'Preferences');
  if strcmp(config.Save_Data_On_Exit, 'yes')
    mifit_File_Save;
  else
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file)), delete(file); end
  end
  mifit_disp([ 'Exiting miFit. Bye bye.' ])
  delete(mifit_fig);
  
function mifit_File_Reset(varargin)
  options.Default     = 'Cancel';
  options.Interpreter = 'tex';
  ButtonName = questdlg({ ...
    '{\fontsize{14}{\color{blue}Reset default configuration ?}}', ...
    'This action will clear the Data set list and the Log file.', ...
    'Selecting "Factory settings" also resets the Preferences in', ...
    prefdir, ...
    '{\bf{Reset now ?}}'}, 'miFit: Reset ?', ...
    'Reset', 'Cancel', 'Factory settings', ...
    options);
  if ~strcmp(ButtonName, 'Cancel')
    if ~strcmp(ButtonName, 'Reset') % Factory settings
      file = fullfile(prefdir, [ mfilename '.ini' ]);
      if ~isempty(dir(file)), delete(file); end
      mifit_Load_Preferences();
      mifit_Apply_Preferences();
    end
    mifit_Edit_Select_All([], 1);
    mifit_Edit_Delete();
    
    file = fullfile(prefdir, [ mfilename '.log' ]);
    if ~isempty(dir(file)), delete(file); end
  end
  
function mifit_File_Log(varargin)
  file = fullfile(prefdir, [ mfilename '.log' ]);
  if ~isempty(dir(file))
    edit(file);
  end
  
% Edit menu ********************************************************************

function mifit_Edit_Undo(varargin)
% set the Data stack to the previous state from History
  mifit_History_pull();
  mifit_List_Data_UpdateStrings();

function mifit_Edit_Cut(varargin)
% get the selected indices in the List, copy these elements to the clipboard
% and delete the elements. Update the History (in Delete).
  mifit_Edit_Copy(varargin{:});
  mifit_Edit_Delete(varargin{:});

function mifit_Edit_Copy(varargin)
% get the selected indices in the List, copy these elements to the clipboard
% we use 'copy' from Y. Lengwiler as it is pure-Matlab, and extends the limited 
% clipboard function.
  d=mifit_List_Data_pull(); % get selected objects
  if isempty(d), return; end
  x = {};
  for index=1:numel(d)
    % we create a cell which contains the Signal and Axes
    if numel(d) == 1, this = d; else this=d(index); end
    for dim = 1:ndims(this)
      x{end+1} = getaxis(this, dim);
    end
    x{end+1} = this{0};
    % add metadata. This also marks the end of the iData object definition for paste.
    x{end+1} = ...
      sprintf('Title="%s";Label="%s";DisplayName="%s";xlabel="%s";ylabel="%s";title="%s";', ...
      this.Title, this.Label, this.DisplayName, xlabel(this), ylabel(this), title(this));
  end
  copy(x);

function mifit_Edit_Paste(varargin)
% append/copy the data sets from the clipboard to the end of the list
% clipboard can be a file name, an iData Tag/ID or ID in the Stack
% for iData ID, make a copy of the objects
% Update the History
  d=paste();
  if iscell(d)
    % we look for numerical cell elements, and split after each non numerical for 
    % new objects
    D = []; this_datax = {}; this_meta = [];
    for index=1:numel(d)
      num = str2num(d{index});
      if ~isempty(num), this_datax{end+1} = num;
      else
        this_meta = str2struct(d{index});
        this = iData(this_datax{:});
        
        this_datax = {};
        if isstruct(this_meta)
          if isfield(this_meta, 'Title'),  this.Title = this_meta.Title; end
          if isfield(this_meta, 'Label'),  this.Label = this_meta.Label; end
          if isfield(this_meta, 'DisplayName'), this.DisplayName = this_meta.DisplayName; end
          if isfield(this_meta, 'xlabel'), xlabel(this, this_meta.xlabel); end
          if isfield(this_meta, 'ylabel'), ylabel(this, this_meta.ylabel); end
          if isfield(this_meta, 'title'),  title(this,  this_meta.title); end
        end
        D = [ D ; this ];
      end
    end % for
    d = D;
  else
    d = iData(d);
  end
  mifit_List_Data_push(d);

function mifit_Edit_Duplicate(varargin)
% copy and paste selected. Update the history (in Paste).
  d=mifit_List_Data_pull();
  mifit_List_Data_push(copyobj(d));
  
function mifit_Edit_Select_All(hObject, select)
% set the List selected values to all ones
% the select argument can be 0 or 1 to deselect/select all
  hObject = mifit_fig('List_Data_Files');
  items   = get(hObject,'String');
  index_selected = get(hObject,'Value');
  if nargin < 2, select=[]; end
  if (numel(index_selected) == numel(items) && isempty(select)) ...
    || (~isempty(select) && ~select), index_selected = [];
  elseif isempty(select) || (~isempty(select) && select), index_selected=1:numel(items); end
  set(hObject,'Value', index_selected);

function mifit_Edit_Delete(varargin)
% delete selected
% Update the History
  fig = mifit_fig;
  hObject        = mifit_fig('List_Data_Files');
  index_selected = get(hObject,'Value');
  Data = getappdata(fig, 'Data');
  if numel(Data) > 1 && numel(index_selected) < numel(Data)
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

function mifit_Data_Plot(varargin)
  d = mifit_List_Data_pull;
  f=figure;
  subplot(d);
  
function mifit_Data_Export(varargin)
  d = mifit_List_Data_pull;
  if ~isempty(d)
    save(d, 'gui');
  end
  
function mifit_Data_Table(varargin)
  d = mifit_List_Data_pull;
  for index=1:numel(d)
    handle = edit(d(index), 'editable');
    set(handle, 'DeleteFcn', @mifit);
  end

function mifit_Data_View(varargin)
  d = mifit_List_Data_pull;
  for index=1:numel(d)  
    if ~isempty(d(index).Source) && ~isdir(d(index).Source)
      try
        edit(d(index).Source)
      end
    end
  end
  
function mifit_Data_Properties(varargin)
  disp('mifit_Data_Properties: should display properties from "disp" and allow to re-assign signal, axes, define new aliases...')
  
function mifit_Data_History(varargin)
  d = mifit_List_Data_pull();
  for index=1:numel(d)
    if ~isempty(d(index).Source) && ~isdir(d(index).Source)
      commandhistory(d(index));
    end
  end
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

function mifit_List_Data_Files(varargin)
% called when clicking on the listbox
disp mifit_List_Data_Files
varargin
%  mifit_Data_Plot();

function mifit_List_Data_push(d)
% put a new data set at the end of the stack
  if isempty(d),       return; end
  if ~isa(d, 'iData'), return; end
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
  
  % Update the History with the new stack
  mifit_History_push;
  mifit_disp('Importing into List:')
  mifit_disp(char(d))
  
function d=mifit_List_Data_pull(varargin)
% get the selected Data List
% return the selected objects
  hObject = mifit_fig('List_Data_Files');
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
  % display message, and log it
  
  if size(message,1) > 1
    disp(message);
  else
    disp([ mfilename ': ' message ]);
  end
  file = fullfile(prefdir, [ mfilename '.log' ]);
  fid = fopen(file, 'a+');
  if fid == -1, return; end
  fprintf(fid, '[%s] %s\n', datestr(now), message);
  fclose(fid);

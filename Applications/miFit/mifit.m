function varargout = mifit(varargin)
% Version: $Date$
% (c) E.Farhi, ILL. License: EUPL.

% notes:
% saveas: display dialogue to select what to export:
% * data sets
% * parameters
% * fit model

% plot: plot model/dataset/parameters ? -> preferences
% will add toolbar instead of uicontrols for fast processing

%  function config = mifit_Load_Preferences
%  function mifit_Save_Preferences(config)
%  function mifit_File_New(handle)
%  function mifit_File_Open(handle)
%  function mifit_File_Save(varargin)
%  function mifit_File_Saveas(varargin)
%  function mifit_File_Print(varargin)
%  function mifit_File_Preferences(varargin)
%  function mifit_Apply_Preferences
%  function mifit_File_Exit(varargin)
%  function mifit_File_Reset(varargin)
%  function mifit_File_Log(varargin)
%  function mifit_Edit_Undo(varargin)
%  function mifit_Edit_Cut(varargin)
%  function mifit_Edit_Copy(varargin)
%  function mifit_Edit_Paste(varargin)
%  function mifit_Edit_Duplicate(varargin)
%  function mifit_Edit_Select_All(hObject, select)
%  function mifit_Edit_Delete(varargin)
%  function mifit_History_pull
%  function mifit_History_push
%  function mifit_Data_Plot(varargin)
%  function mifit_Data_Export(varargin)
%  function mifit_Data_Table(varargin)
%  function mifit_Data_View(varargin)
%  function mifit_Data_Properties(varargin)
%  function mifit_Data_History(varargin)
%  function mifit_Tools_About(fig)
%  function mifit_List_Data_Files(varargin)
%  function mifit_List_Data_push(d)
%  function mifit_List_Data_pull(varargin)
%  function mifit_List_Data_UpdateStrings
%  function mifit_disp(message)

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
          if strcmpi(action,'identify'), varargout{1} = []; return; end
          if any(strcmpi(action,{'pull','data'}))
            out = getappdata(fig, 'Data');
          else
            try
              feval([ 'mifit_' action ], varargin{2:end});
            catch ME
              mifit_disp([ '[mifit] Unknown action "' action '"' ])
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
        elseif isa(varargin{1}, 'iData')
          mifit_List_Data_push(varargin{1});
        elseif isa(varargin{1}, 'iFunc') || isa(varargin{1}, 'sw') || isa(varargin{1}, 'spinw')
          % TODO: would push new Model
          % mifit_List_Data_push(varargin{1});
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

% Preferences I/O --------------------------------------------------------------




% -------------------------------------------------------------------------
% Callbacks
% -------------------------------------------------------------------------
% usual callback arguments are : Callback(gcbo,eventdata)

% File menu ********************************************************************
function mifit_File_New(handle)
% File/New menu item
  d = iData(zeros(5)); % create an empty Data set;
  mifit_disp([ '[File_New] Editing an empty Data set. ' ...
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
  
  file = fullfile(prefdir, [ mfilename '.mat' ]);
  mifit_File_Saveas(file);
  
function mifit_File_Saveas(varargin)
% save the application configuration into specified file
% Data sets and Model parameters into a mifit.mat file
  fig        = mifit_fig;
  Data       = getappdata(fig, 'Data');
  Models     = getappdata(fig, 'Models');
  Optimizers = getappdata(fig, 'Optimizers');
 
  if nargin == 1 && ischar(varargin{1})
    file = varargin{1};
  else
    filterspec = { '*.mat','MAT-files (*.mat)'};
    [filename, pathname] = uiputfile(filterspec, 'Save All miFit Data sets and models as', [ mfilename '.mat' ]);
    if isequal(filename,0) || isequal(pathname,0)
      return
    end
    file = fullfile(pathname, filename);
  end
  mifit_disp([ '[File_Saveas] Saving Data sets/Models into ' file ]);
  builtin('save', file, 'Data','Models','Optimizers');
  
function mifit_File_Print(varargin)
% print the interface. Not very nice. can we think of something better ?
% perhaps we can generate an HTML report in Saveas HTML ?
  fig = mifit_fig;
  printdlg(fig);
  % alternative: File_Saveas_HTML in tmpfile, then open that file for printing.
  disp([ mfilename ': File_Print: TODO: save all Data sets as HTML with model and parameters' ])
  disp('then open it with web for printing');
  
function mifit_File_Preferences(varargin)
% open Preferences dialogue
% set directories to search for Models
% set FontSize (and update all Fonts in figure)
% set Save on exit
% save Preferences on dialogue close
  fig = mifit_fig;
  config = getappdata(mifit_fig, 'Preferences');
  prompt = {'Font size [10-36]','Save Data sets on Exit [yes/no]','Store Models when creation time is longer than [sec, 0:always, Inf:never, default=3]'};
  if ~isfield(config, 'FontSize'),          config.FontSize=12; end
  if ~isfield(config, 'Save_Data_On_Exit'), config.Save_Data_On_Exit='yes'; end
  if ~isfield(config, 'Store_Models'),      config.Store_Models=3; end
  defaultanswer = { num2str(config.FontSize), config.Save_Data_On_Exit, num2str(config.Store_Models) };
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
  answer{3} = str2double(answer{3});
  if isfinite(answer{3}) 
    config.Store_Models=answer{3}; end
  setappdata(mifit_fig, 'Preferences', config);
  mifit_Apply_Preferences;
  mifit_Save_Preferences(config);

function mifit_File_Exit(varargin)
% Quit and Save Data
  config = getappdata(mifit_fig, 'Preferences');
  if isfield(config, 'Save_Data_On_Exit') && strcmp(config.Save_Data_On_Exit, 'yes')
    mifit_File_Save;
  else
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file)), delete(file); end
  end
  mifit_disp([ '[Exit] Exiting miFit. Bye bye.' ])
  delete(mifit_fig);
  
function mifit_File_Reset(varargin)
  options.Default     = 'Cancel';
  options.Interpreter = 'tex';
  ButtonName = questdlg({ ...
    '{\fontsize{14}{\color{blue}Reset default configuration ?}}', ...
    'Selecting "Reset" clears the Data set list and the Log file.', ...
    'Selecting "Factory settings" also resets the Models and Preferences in', ...
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
      setappdata(mifit_fig, 'Models',[]);
      setappdata(mifit_fig, 'Optimizers',[]);
      file = fullfile(prefdir, [ mfilename '.mat' ]);
      if ~isempty(dir(file)), delete(file); end
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
  
% Data menu ********************************************************************

function mifit_Data_Plot(varargin)
  d = mifit_List_Data_pull;
  f=figure;
  subplot(d);
  
function mifit_Data_Fit(varargin)
  d = mifit_List_Data_pull;
  % TODO: it is desirable to handle constraints (min/max/fix)
  p=fits(d, '', '', 'OutputFcn=fminplot;Display=iter');  % with assigned models or gaussians
  
function mifit_Data_Saveas(varargin)
  d = mifit_List_Data_pull;
  if ~isempty(d)
    save(d, 'gui');
  end
  
function mifit_Data_Table(varargin)
  d = mifit_List_Data_pull;
  config = getappdata(mifit_fig, 'Preferences');
  
  for index=1:numel(d)
    handle = edit(d(index), 'editable');
    set(handle, 'FontSize', config.FontSize);
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
% TODO
  disp([ mfilename ': Data_Properties: TODO: should display properties from "disp" and allow to re-assign signal, axes, define new aliases...' ])
% Data_Properties ?
% re-assign data set signal, axes, ... to aliases/new ones
% display statistics
  
function mifit_Data_History(varargin)
  d = mifit_List_Data_pull();
  for index=1:numel(d)
    [c,fig]=commandhistory(d(index));
  end
  
function mifit_Data_Math_Unary(varargin)
  % TODO
  disp([ mfilename ': Data_Math_Unary: TODO' ])
  
function mifit_Data_Math_Binary(varargin)
  % TODO
  disp([ mfilename ': Data_Math_Binary: TODO' ])
  
function mifit_Data_Math(varargin)
  % TODO
  disp([ mfilename ': Data_Math (others): TODO' ])

% Models and Optimizers menu ***************************************************

function mifit_Models_Add(varargin)
  % TODO
  disp([ mfilename ': Models_Add: TODO' ])
  % * Add from file... (JSON, M, YAML, MAT)

function mifit_Models_Edit(varargin)
  % TODO
  disp([ mfilename ': Models_Edit: TODO' ])
  % create new Models after edition
  
function mifit_Models_Plot(varargin)
  % TODO
  disp([ mfilename ': Models_Plot: TODO' ])
  
function mifit_Models_Plot_Parameters(varargin)

function mifit_Models_Export_Parameters(varargin)
  % * Export to file... (JSON, M, YAML, MAT...)

function mifit_Models_View_Parameters(varargin)
  % get 1st selected Model from Data set or Models menu current choice
  % Display a uitable with columns:
  % [ Parameters | ParameterValues | constraints.fixed | constraints.min | constraints.max ]
  
function mifit_Models_Add_Expression(varargin)
    % * 4D TAS convolution        -> in Models Transformation/operations
  % * Powder average 4D -> 2D   -> in Models Transformation/operations

% set optimizer configuration -> contextual dialogue in Model_Parameters uitable ?
  
% Tools menu *******************************************************************
  
function mifit_Tools_Help_Loaders(varargin)
  doc(iData,'Loaders');
  
function mifit_Tools_Help_Models(varargin)
  doc(iData,'Models');
  
function mifit_Tools_Help_Optimizers(varargin)
  doc(iData,'Optimizers');
  
% List Data and Stack management ***********************************************

function mifit_List_Data_Files(varargin)
% called when clicking on the listbox

%  mifit_Data_Plot();
  
function mifit_List_Data_UpdateStrings
  % update the List labels
  fig = mifit_fig;
  
  Data = getappdata(fig, 'Data');
  
  hObject        = mifit_fig('List_Data_Files');
  list           = {};
  if numel(Data) > 1
    for index=1:numel(Data)
        list{end+1} = char(Data(index));
    end
  else
    list{end+1} = char(Data);
  end
  set(hObject,'String', list, 'Value', []);


  

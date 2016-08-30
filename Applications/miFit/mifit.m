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
    mifit_disp('[Init] Welcome to miFit ! ********************************************')
    fig = openfig(mfilename);
    
    % get Models/Optimizers menu handles
    hmodels = mifit_fig('Menu_Model');
    hoptim  = mifit_fig('Menu_Optimizers');
    
    % load Preferences
    mifit_Load_Preferences;
    mifit_Apply_Preferences;
    % create the AppData default values
    setappdata(fig, 'Data',    []);
    setappdata(fig, 'History', {});
    setappdata(fig, 'Models',  {});
    
    % Display welcome dialog during menu build
    h = mifit_Tools_About(fig);
    contrib = textwrap({ version(iData,2) },80);
    for index=1:numel(contrib)
      mifit_disp([ '  ' contrib{index} ])
    end
    
    % get the list of Models and Optimizers
    models = [];
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file))
      try
        d = load(file);
        if isfield(d, 'Models')
          mifit_disp([ '[Init] Loading Model library from ' file ]);
          models = d.Models;  % contains the callback and the label members
        end
      end
    end
    if isempty(models)
      mifit_disp([ '[Init] Building Optimizer and Model library from ' file '. Be patient (only once)...' ]);
      [optimizers,models,filenames] = fits(iFunc);
      % build Models.callback and Models.label. We will instantiate at callback.
      models = struct('callback',filenames,'label',get(models, 'Name'));
    else
      optimizers = fits(iFunc);
    end
    
    % fill Models menu
    mifit_Models_Add_Entry(models);
    
    % fill Optimizers menu
    if ~isempty(optimizers) && iscell(optimizers)
        mifit_disp([ '[Init] Initializing ' num2str(numel(optimizers)) ' Optimizers ...' ]);
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
              % TODO: must add callback to assign optimizer
              uimenu(hoptim, 'Label', algorithm, 'UserData', f{1});
            end
        end
    end
    
    % create the AppData Data Stack
    setappdata(fig, 'Models',    models);
    setappdata(fig, 'Optimizers',optimizers);
    
    % Load the previous Data sets containing Model Parameters (when a fit was performed)
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file))
      try
        d = load(file);
        if isfield(d, 'Data')
          mifit_disp([ '[Init] Loading Data sets from ' file ]);
          mifit_List_Data_push(d.Data);
        end
      end
    end
    file = fullfile(prefdir, [ mfilename '.log' ]);
    mifit_disp([ '[Init] Log file is ' file ]);

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
      mifit_disp([ '[Load_Preferences] Loading Preferences from ' file ]);
    end
  end
  if isempty(config)
    % default configuration
    config.FontSize         = max(12, get(0,'defaultUicontrolFontSize'));
    config.Save_Data_On_Exit= 'yes';
    config.Store_Models     = 10;  % time required for creation. Store when > 0:always, Inf=never
  end
  setappdata(mifit_fig, 'Preferences', config);
  set(0,'defaultUicontrolFontSize', config.FontSize);

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
    mifit_disp([ '[Save_Preferences] Saving Preferences into ' filename ]);
  end

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
  if isempty(Data), return; end
 
  if nargin == 1 && ischar(varargin{1})
    file = varargin{1};
  else
    filterspec = { '*.mat','MAT-files (*.mat)'};
    [filename, pathname] = uiputfile(filterspec, 'Save All miFit Data sets as', [ mfilename '.mat' ]);
    if isequal(filename,0) || isequal(pathname,0)
      return
    end
    file = fullfile(pathname, filename);
  end
  mifit_disp([ '[File_Saveas] Saving Data sets into ' file ]);
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
  prompt = {'Font size [10-36]','Save Data sets on Exit [yes/no]','Store Models when creation time is longer than [sec, 0:always, Inf:never, default=10]'};
  if ~isfield(config, 'FontSize'),          config.FontSize=12; end
  if ~isfield(config, 'Save_Data_On_Exit'), config.Save_Data_On_Exit='yes'; end
  if ~isfield(config, 'Store_Models'),      config.Store_Models=10; end
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
  
function mifit_Apply_Preferences
  fig = mifit_fig;
  config = getappdata(fig, 'Preferences');
  % change all Font Size
  h=[ findobj(fig, 'Type','uicontrol') ; ...
      findobj(fig, 'Type','axes') ; findobj(fig, 'Type','text') ; ...
      findobj(fig, 'Type','uipanel') findobj(fig, 'Type','uitable') ];
  set(h, 'FontSize', config.FontSize);
  % for uimenu, could we use tip given by Y Altman 
  % <https://fr.mathworks.com/matlabcentral/newsreader/view_thread/148095>
  % h = findobj(fig, 'Type','uimenu');

function mifit_File_Exit(varargin)
% Quit and Save Data
  config = getappdata(mifit_fig, 'Preferences');
  if strcmp(config.Save_Data_On_Exit, 'yes')
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
  if ~isempty(Data)
    History{end+1} = Data;
    if numel(History) > 10, History(1:(end-9)) = []; end
    setappdata(fig, 'History', History);
  end
  
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

function mifit_Data_AssignModel(varargin)
  model = get(varargin{1},'UserData'); % an iFunc or char/cellstr stored into UserData of the menu item.
  
  if iscellstr(model), model = char(model)'; end
  if ischar(model)
    try
      % create the Model from the stored expression
      tstart   = tic;
      model    = eval(model);
      telapsed = toc(tstart);
      if telapsed > 10
        % TODO: push Model into the Models menu (to the Deck/Library/Stack)
      end
    catch
      mifit_disp([ '[Data_AssignModel] Invalid Model expression ' model '. Skipping.' ]);
      return
    end
  elseif isa(model, 'iFunc')
    model = copyobj(model); % to get a new ID
  end
  
  setappdata(mifit_fig, 'CurrentModel', model);  % store current selected Model
  % get selected Data sets indices in List
  index_selected = get(mifit_fig('List_Data_Files'),'Value');
  D = getappdata(mifit_fig, 'Data');  % all data sets
  if numel(D) == 0 || isempty(index_selected), 
    mifit_disp([ 'Selected Model "' model.Name '".' ]);
    figure; plot(model);
    return; 
  end
  
  mifit_disp([ 'Assigning Model "' model.Name '" to ' num2str(numel(index_selected)) ' Data set(s).' ]);
  mifit_History_push();
  if numel(D) > 1
    for index=index_selected(:)'
      D(index) = setalias(D(index), 'Model', model);
      mifit_disp(char(D(index)));
    end
  else
    D = setalias(D, 'Model', model);
    mifit_disp(char(D));
  end
  setappdata(mifit_fig, 'Data', D);
  
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
 
  
function mifit_Models_Add_Entry(model)
  % add a new Model (iFunc or cellstr) into the Models sub-menus.
  % the input argument should be a structure with members:
  %   model.callback
  %   model.label
  % or an iFunc array
  % or a cellstr (expressions to evaluate at callback)
  
  % usage:
  % mifit_Models_Add_iFunc(gauss)         - must add a Gaussian entry in sub-menu, with callback to 'copyobj(gauss)'
  % mifit_Models_Add_iFunc('gauss+lorz')  - must add a 'gauss+lorz' entry in sub-menu, with callback to 'gauss+lorz'
  
  if nargin == 0 || all(isempty(model)), return; end
  
  % check different type of input arguments Models
  if ischar(model), model = cellstr(model); end
  
  % handle array of entries
  if numel(model) > 1  % model is a cell or array of entries
    for index=1:numel(model)
      if iscell(model), mifit_Models_Add_Entry(model{index});
      else mifit_Models_Add_Entry(model(index));
      end
    end
    return
    
  % model is a single struct with 'callback' as an array
  elseif isstruct(model) && isfield(model,'callback') && ~ischar(model.callback) && numel(model.callback) > 1
    for index=1:numel(model.callback)
      % we create a model struct with only one (callback,label) item
      if isfield(model,'label') && iscell(model.label) && numel(model.label) == numel(model.callback), 
        this.label=model.label{index}; 
      else 
        this.label=''; 
      end
      if iscell(model.callback), this.callback = model.callback{index};
      else                       this.callback = model.callback(index); end
      mifit_Models_Add_Entry(this);
    end
    return
  end
  
  % handle single entry
  % determine
  %   label: a string to display for the menu entry
  %   callback: a char/cellstr (expression) or iFunc model, stored in the menu 
  %     item UserData for assignement when selected
  %   dim:   dimensionality (to identify the sub-menu)
  
  label = ''; dim = [];
  if all(isempty(model)), return; end
  if iscellstr(model)
    callback = char(model)';
  elseif isstruct(model) && isfield(model,'callback')
    callback = model.callback;
    if isfield(model, 'label'), label = model.label; end
  elseif isa(model, 'iFunc')
    callback = model;
    dim      = callback.Dimension;
    label    = [ callback.Name ' [' callback.Tag ']' ];
  elseif ischar(model)
    callback = model;
  else
    mifit_disp([ '[Models_Add_Entry] Invalid Model type ' class(model) '. Must be iFunc/cell/char/struct with "callback" member. Skipping.' ]);
  end
  
  if isempty(dim)  % this is a model creator (new instance) from expression
    % we assume this is an expression and try to evaluate it as such
    try
      modelF    = feval(callback, 'identify');
    catch
      mifit_disp([ '[Models_Add_Entry] Invalid Model expression ' callback '. Skipping.' ]);
      return
    end
    dim   = modelF.Dimension;
    if isempty(label), label = [ '"' callback '" = ' modelF.Name ]; end
  end

  % determine the name of the sub-menu to use
  if isempty(dim), dim = -1; end  % will use Others
  if dim > 0, model_submenu_name = sprintf('%dD', dim);
  else        model_submenu_name = 'Others'; end
  
  % create the sub-menu, if missing
  hmodels = mifit_fig('Menu_Model');
  submenu_handle = findobj(hmodels, 'Type','uimenu', 'Label', model_submenu_name);
  if isempty(submenu_handle)
    submenu_handle = uimenu(hmodels, 'Label', model_submenu_name);
  end
  
  % check if the model entry already exists
  children = findobj(submenu_handle, 'Label', label);
  if isempty(children)
    children = findobj(submenu_handle, 'UserData', callback);
  end
  if ~isempty(children)
    mifit_disp([ '[Models_Add_Entry] ' label ' is already in the list of usable Models. Skipping.' ]);
    return; 
  end
    
  uimenu(submenu_handle, 'Label', label, 'UserData', callback, ...
                'CallBack', 'mifit(''Data_AssignModel'',gcbo)');
                
  % store the entry in the appdata
  % Models is a cell of entries which can be:
  %   struct.callback
  %   struct.callback and struct.label
  %   iFunc
  %   char (expression)
  Models = getappdata(mifit_fig, 'Models');
  Models{end+1} = struct('callback',callback, 'label',label);
  setappdata(mifit_fig, 'Models',Models);

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
  g=findobj(h, 'type','uicontrol');
  config = getappdata(mifit_fig, 'Preferences');
  set(g,'fontsize',config.FontSize);
  if ~isempty(fig)
    setappdata(fig, 'handle_About', h);
  end
  
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
  
function [d, index_selected]=mifit_List_Data_pull(varargin)
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
  if numel(Data) > 1
    for index=1:numel(Data)
        list{end+1} = char(Data(index));
    end
  else
    list{end+1} = char(Data);
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
  

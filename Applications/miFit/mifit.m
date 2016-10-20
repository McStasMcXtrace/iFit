function varargout = mifit(varargin)
% miFit: a user interface to iFit
%
% data = mifit('data')  retrieves data sets from the interface
% mifit('filename')     imports the file into a new Data set/Model
% mifit(iData_object)   add the iData object into the interface Data stack
% mifit(iFunc_object)   add the iFunc Model into the interface Models menu
%
% Version: $Date$
% (c) E.Farhi, ILL. License: EUPL.

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
%  function mifit_Data_Fit(varargin)
%  function mifit_Help_About(fig)
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
    else set(fig, 'NextPlot','new','HandleVisibility','callback');
    end
    out = fig;

    if ~isempty(varargin)
        if nargin == 2 && isa(varargin{1}, 'dndcontrol') && isstruct(varargin{2})
          % called as a callback/event from Drag-n-Drop
          evt = varargin{2};
          switch evt.DropType
          case 'file'
              for n = 1:numel(evt.Data)
                  mifit(deblank(evt.Data{n}));
              end
          case 'string'
              lines = textscan(evt.Data, '%s','Delimiter','\n');
              mifit(deblank(lines{1}));
          end
        elseif ischar(varargin{1}) && isempty(dir(varargin{1})) ...
          && ~any(strncmp(varargin{1}, {'http:','ftp:/','file:','https'}, 5)) % a function/action to call ?
          % callback with varargin{1} == 'action'
          action = varargin{1};
          if strcmpi(action,'identify'), varargout{1} = []; return; end
          if any(strcmpi(action,{'pull','data'}))
            out = getappdata(fig, 'Data');
          else
            try
              feval([ 'mifit_' action ], varargin{2:end});
            catch ME
              mifit_disp([ '[mifit] Error for action "' action '"' ])
              rethrow(ME)
            end
          end
        elseif ishandle(varargin{1})    % a CallBack function
          h = varargin{1}; d = [];
          UserData = get(h,'UserData');
          if isfield(UserData,'modified') && ~UserData.modified
            d = [];
          elseif strcmpi(get(h,'Type'),'uitable')
            % get data from a UITable
            try
                d = edit(iData, h);
            catch
                d = iData(get(h, 'Data'));
                d.Title = get(h, 'Tag');
            end
            % now push 'd' into the Stack, except for 'empty'/non modified Table
          else
            d = iData(h);
          end
          if ~isempty(d)
            mifit_List_Data_push(d);
          end
        elseif isa(varargin{1}, 'iData')
          mifit_List_Data_push(varargin{1});
        elseif isa(varargin{1}, 'iFunc') || isa(varargin{1}, 'sw') || isa(varargin{1}, 'spinw')
          mifit_Models_Add_Entry(iFunc(varargin{1}));
        else  % import files/other stuff
          d = [];
          if ischar(varargin{1}) && ~isempty(dir(varargin{1}))
            % test if we import a saved MAT file
            if ~isdir(varargin{1})
              file = varargin{1};
              try
                d = load(file);
                if isfield(d, 'Data'), d=d.Data; else d=[]; end
              catch
                d = [];
              end
              % test if we import a config.ini file
              if isempty(d)
                try
                  config = mifit_Load_Preferences(file);
                catch
                  config = [];
                end
                if ~isempty(config)
                  mifit_Apply_Preferences(config);
                  return
                end
              end
            end
          end
          if isempty(d)
            d = iData(varargin{:});
          end
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
% File/New menu item: open a table to edit a new data set
  d = iData(zeros(5)); % create an empty Data set;
  mifit_disp([ '[File_New] Editing an empty Data set. ' ...
    'Close the window to retrieve its content as a new Data set into miFit. ' ...
    'Use the Contextual menu for Copy/Paste/Resize.' ]);
  handle = edit(d, 'editable');
  % set a DeletedFcn so that the content can be retrieved into miFit when
  % closing.
  set(handle, 'DeleteFcn', @mifit);
  config = getappdata(mifit_fig, 'Preferences');
  set(handle, 'FontSize', config.FontSize);
  
function mifit_File_Open(handle)
% File/Open menu item: shown fileselector and open file
  d = iData('');  % open file selector, and import files
  % push that data onto the List
  mifit_List_Data_push(d);
  
function mifit_File_Save(varargin)
% File/Save: save Data sets and Model parameters into a mifit.mat file
  
  file = fullfile(prefdir, [ mfilename '.mat' ]);
  mifit_File_Saveas(file);
  
function mifit_File_Saveas(varargin)
% File/saveas: save as a MAT with user filename
% save the application configuration into specified file
% Data sets and Model parameters into a mifit.mat file
  fig        = mifit_fig;
  Data       = getappdata(fig, 'Data');
  Models     = getappdata(fig, 'Models');
  Optimizers = getappdata(fig, 'Optimizers');
  CurrentOptimizer         = getappdata(fig, 'CurrentOptimizer');
  CurrentOptimizerConfig   = getappdata(fig, 'CurrentOptimizerConfig');
  CurrentOptimizerCriteria = getappdata(fig, 'CurrentOptimizerCriteria');
  CurrentModel             = getappdata(fig, 'CurrentModel');
 
  if nargin == 1 && ischar(varargin{1})
    file = varargin{1};
  else
    filterspec = { '*.mat','MAT-files (*.mat)'};
    [filename, pathname] = uiputfile(filterspec, 'Save All miFit Data sets and Models as', [ mfilename '.mat' ]);
    if isequal(filename,0) || isequal(pathname,0)
      return
    end
    file = fullfile(pathname, filename);
  end
  builtin('save', file, 'Data','Models','Optimizers','CurrentOptimizer', ...
    'CurrentOptimizerConfig', 'CurrentOptimizerCriteria');
  d = dir(file);
  mifit_disp([ '[File_Saveas] Saved Data sets/Models into ' file ' [size: ' num2str(round(d.bytes/1024)) ' kb]' ]);
  
function mifit_File_Print(varargin)
% File/Print: print the interface. 
% generate an HTML report and display in browser for printing.
  d=mifit_List_Data_pull(); % get selected objects
  if all(isempty(d)), return; end
  filename = [ tempname '.html' ];
  mifit_disp([ '[File_Print] Exporting Data sets to HTML ' filename ' for printing...' ]);
  save(d, filename, 'html data');
  webbrowser(filename,'system');  % tries to open with the system browser
  
function mifit_File_Preferences(varargin)
% File/Preferences: open Preferences dialogue
% set directories to search for Models
% set FontSize (and update all Fonts in figure)
% set Save on exit
% save Preferences on dialogue close
  fig = mifit_fig;
  config = getappdata(mifit_fig, 'Preferences');
  if ~isfield(config, 'FontSize'),          config.FontSize=12; end
  if ~isfield(config, 'Save_Data_On_Exit'), config.Save_Data_On_Exit='yes'; end
  if ~isfield(config, 'Store_Models'),      config.Store_Models=3; end
  if ~isfield(config, 'History_Level'),     config.History_Level=10; end
  if ~isfield(config, 'Fit_Verbose'),       config.Fit_Verbose='yes'; end
  options.Name       = [ mfilename ': Preferences' ];
  options.ListString = {'FontSize Font size [10-36]', ...
    'Save_Data_On_Exit Save Data sets on Exit [yes/no]', ...
    'Store_Models Store Models when creation time is longer than [sec, 0:always, Inf:never, default=3]', ...
    'History_Level Undo levels to keep [2-50, reduce if you handle large/many data sets]', ...
    'Fit_Verbose Verbosity when Fitting [yes shows criteria, parameters and final distributions]' };
  options.FontSize   = config.FontSize;
  options.TooltipString = sprintf([ 'Modify the %s Preferences.\n' ...
      'You can specify additional menu items by entering a cell (pairs)\n' ...
      '* Menu_<Label> = {''Item_Label'',''Command'', ...}' ], mfilename);
  config1 = structdlg(config, options);
  if isempty(config1), return; end
  config = config1;
  % set new Preferences
  config.FontSize=min(max(config.FontSize, 10),36);
  config.History_Level=min(max(config.History_Level, 2),50);
  setappdata(mifit_fig, 'Preferences', config);
  mifit_Apply_Preferences;
  mifit_Save_Preferences(config);

function mifit_File_Exit(varargin)
% File/Exit: Quit and Save Data
  config = getappdata(mifit_fig, 'Preferences');
  if isfield(config, 'Save_Data_On_Exit') && strcmp(config.Save_Data_On_Exit, 'yes')
    mifit_File_Save;
  else
    file = fullfile(prefdir, [ mfilename '.mat' ]);
    if ~isempty(dir(file)), delete(file); end
  end
  mifit_disp([ '[Exit] Exiting miFit. Bye bye.' ]);
  delete(mifit_fig('mifit_View_Parameters'));
  delete(mifit_fig);
  
function mifit_File_Reset(varargin)
% File/Reset: clear data sets, user models and Preferences
  options.Default     = 'Cancel';
  options.Interpreter = 'tex';
  ButtonName = questdlg({ ...
    '{\fontsize{14}{\color{blue}Reset default configuration ?}}', ...
    'Selecting "Reset" clears the Data set list, History and the Log file.', ...
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
    mifit_Edit_Select_All([], true);
    mifit_Edit_Delete();
    setappdata(fig, 'History', []);

    file = fullfile(prefdir, [ mfilename '.log' ]);
    if ~isempty(dir(file)), delete(file); end
  end
  
function mifit_File_Log(varargin)
% File/Log: show mifit log file
  file = fullfile(prefdir, [ mfilename '.log' ]);
  if ~isempty(dir(file))
    edit(file);
  end
  
% Edit menu ********************************************************************

function mifit_Edit_Undo(varargin)
% Edit/Undo: set the Data stack to the previous state from History
  mifit_History_pull();
  mifit_List_Data_UpdateStrings();

function mifit_Edit_Cut(varargin)
% Edit/Cut: get the selected indices in the List, copy these elements to the clipboard
% and delete the elements. Update the History (in Delete).
  mifit_Edit_Copy(varargin{:});
  mifit_Edit_Delete(varargin{:});

function mifit_Edit_Copy(varargin)
% Edit/Copy: get the selected indices in the List, copy these elements to the clipboard
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
% Edit/Paste: append/copy the data sets from the clipboard to the end of the list
% clipboard can be a file name, an iData Tag/ID or ID in the Stack
% for iData ID, make a copy of the objects
% Update the History
  d=paste();
  if iscell(d)
    % first look for file names and URL's
    D = [];
    for index=1:numel(d)
      if ~isempty(dir(d{index})) || any(strncmp(d{index}, {'http:','ftp:/','file:','https'},5))
        D = [ D ; iData(d{index}) ];
        d{index} = [];
      end
    end
    if ~isempty(D), mifit_List_Data_push(D); end
    % we look for numerical cell elements, and split after each non numerical for 
    % new objects
    D = []; this_datax = {}; this_meta = [];
    for index=1:numel(d)
      if isempty(d{index}), continue; end
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
% Edit/Duplicate: copy and paste selected. Update the history (in Paste).
  d=mifit_List_Data_pull();
  mifit_List_Data_push(copyobj(d));

function mifit_Edit_Delete(varargin)
% Edit/Delete: delete selected
% Update the History
  fig = mifit_fig;
  hObject        = mifit_fig('List_Data_Files');
  index_selected = get(hObject,'Value');
  if isempty(index_selected), return; end
  Data = getappdata(fig, 'Data');
  if numel(Data) == 0, return; end
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
% Data/Plot
  if nargin && isa(varargin{1}, 'iData'), d=varargin{1}; 
  else d = mifit_List_Data_pull; end
  if all(isempty(d)), return; end
  f=figure;
  subplot(d,'light transparent grid tight');
  
function mifit_Data_Saveas(varargin)
% Data/Saveas: export selected data sets
  d = mifit_List_Data_pull;
  if numel(d)
    save(d, 'gui');
  end
  
function mifit_Data_Table(varargin)
% Data/view Table: view data sets as tables
  d = mifit_List_Data_pull;
  config = getappdata(mifit_fig, 'Preferences');
  
  for index=1:numel(d)
    handle = edit(d(index), 'editable');
    set(handle, 'FontSize', config.FontSize);
    set(handle, 'DeleteFcn', @mifit);
  end

function mifit_Data_View(varargin)
% Data/view Source: edit source file when exists
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

function mifit_Data_Search(varargin)
% Data/Search: search all Data sets for a token
  prompt = {'Enter the item to search for (string, case sensitive):'};
  name   = 'miFit: Search in all data sets';
  numlines=1;
  answer =inputdlg(prompt,name,numlines);
  if isempty(answer),    return; end  % cancel
  if isempty(answer{1}), return; end  % empty pattern
  
  D=getappdata(mifit_fig, 'Data');
  
  mifit_disp([ '[Data search] Data sets mentioning "' answer{1} '" ...' ]);
  % search for token in strings
  match = strfind(D, answer{1}, 'case');
  index_selected = ~cellfun(@isempty, match);
  % search for token in field names
  match = findfield(D, answer{1}, 'case');
  index_selected = index_selected | ~cellfun(@isempty, match);
  index_selected = find(index_selected);
  
  if isempty(index_selected)
    mifit_disp('None');
  else
    mifit_disp(char(D(index_selected)));
    mifit_Edit_Select_All([], index_selected);
  end
  
  
function mifit_Data_History(varargin)
% Data/History: show Data sets command history
  d = mifit_List_Data_pull();
  for index=1:numel(d)
    [c,fig]=commandhistory(d(index));
  end
  
function mifit_Data_Math_Unary(varargin)
% Data/Transform unary: apply a unary operator on Datasets
  
  op = {'abs','acosh','acos','asinh','asin','atanh','atan','ceil','conj','cosh','cos', ...
    'ctranspose','exp','fliplr','flipud','floor','full','imag','isfinite', ...
    'isfloat','isinf','isinteger','islogical','isnan','isnumeric','isreal','isscalar', ...
    'issparse','log10','log','norm','not','permute','real','round','sign','sinh','sin', ...
    'sparse','sqrt','tanh','tan','transpose','uminus','uplus', ...
    'min','mean','median','sort','squeeze',...
    'double','logical','single','sqr'};
    
  mifit_Data_Operator(op);
  
function mifit_Data_Math_Binary(varargin)
% Data/Transform binary: apply a binary operator on Datasets
  
  op = {'combine','conv','convn','eq','ge','gt','le','lt','minus', ...
  'mrdivide','ne','plus','power','rdivide','times','xcorr','interp'};
  
  mifit_Data_Operator(op);
  
function mifit_Data_Math(varargin)
% Data/Transform binary: apply a binary operator on Datasets

  % should contain e.g. conv4d and sqw_powder...
  op = {'cwt','fft','ifft','del2','gradient','diff',...
    'cumtrapz','sum','prod','trapz','cumsum','cumprod',...
    'smooth','kmeans','pca','interp','std','fill','hist'};
  mifit_Data_Operator(op);

function mifit_Data_Operator(op)
% Data/Transform: apply a binary operator on Datasets

  % we display a dialogue to select which operator to apply on the selected Datasets
  [d, index_selected] = mifit_List_Data_pull;
  if isempty(index_selected) || numel(d) == 0, return; end
  
  op = sort(op); % alpha order sort
    
  % show operator selection dialogue
  PromptString = cellstr(char(d));
  if numel(PromptString) > 3
    PromptString=[ PromptString(1:3) ; ' ...' ];
  end
  [selection, ok] = listdlg('ListString', op, 'ListSize',[500 160], ...
    'Name',[ 'miFit: Select operator to apply on ' ...
      num2str(numel(index_selected)) ' data set(s)' ], ...
      'PromptString', PromptString, 'SelectionMode','single');
  if isempty(selection) || ok ~= 1, return; end
  op = op{selection};
 
  % apply operator on Datasets
  mifit_disp([ 'Applying operator "' op '" on data set(s):' ]);
  mifit_disp(char(d))
  try
    d = feval(op, d);
  catch ME
    mifit_disp([ '[Data operation] Failed applying operator "' op '"' ]);
    mifit_disp(getReport(ME));
    return
  end
    
  % upload new Data sets
  mifit_List_Data_push(d);
  
% Models menu ******************************************************************

function [ifuncs, labels,indices,handles] = mifit_Models_GetList(varargin)
% [internal] mifit_Models_GetList: get the list of iFunc models (not expressions) in the Models menu
% indices is the index of 'static' Models in the whole list.
  models = getappdata(mifit_fig,'Models');
  ifuncs = []; labels = {}; indices = []; handles = [];
  for index=1:numel(models)
    this = models{index};
    if ~isempty(this.callback) && isa(this.callback, 'iFunc')
      ifuncs = [ ifuncs this.callback ];
      handles= [ handles this.handle ];
      labels{end+1} = this.label;
      indices(end+1) = index;
    end
  end

function mifit_Models_Load(varargin)
% Models/Import: load new Models
% we pop-up a file selector, read as iFunc, and store it in miFit.
  models = load(iFunc, '');
  if isempty(models), return; end
  mifit_disp([ '[Models] Loading ' num2str(numel(models)) ' new models:' ]);
  mifit_disp(display(models));
  mifit(models);
  
function mifit_Models_Export(varargin)
% Models/Export: save user Models
  % get the list of 'static' iFunc models (which have been created and stored in the Models menu)
  [ifuncs, labels,indices] = mifit_Models_GetList();
  if isempty(indices), return; end
  % pop-up a dialogue box to select those to export, with select all button
  [selection, ok] = listdlg('ListString', labels, 'SelectionMode', 'multiple', ...
    'Name','miFit: Select Models to Export', 'ListSize',[400 160]);
  if isempty(selection) || ok ~= 1, return; end
  % pop-up the iFunc.save export dialogue
  save(ifuncs(selection),'gui');
  
function mifit_Models_Remove(varargin)
% Models/Remove: clear user models
  % select 'static' iFunc models to remove from the menu
   
  % get the list of 'static' iFunc models (which have been created and stored in the Models menu)
  [ifuncs, labels, indices,handles] = mifit_Models_GetList();
  if isempty(indices), return; end
  % pop-up a dialogue box to select those to remove, with select all button
  [selection, ok] = listdlg('ListString', labels, 'SelectionMode', 'multiple', ...
    'Name','miFit: Select Models to Remove', 'ListSize',[400 160]);
  if isempty(selection) || ok ~= 1, return; end
  delete(handles(selection))
  models = getappdata(mifit_fig,'Models');
  models(indices(selection)) = [];
  setappdata(mifit_fig,'Models',models);

function mifit_Models_Edit(varargin)
  % TODO: should be a free input dialogue where we can enter expressions, 
  % and insert existing objects reference
  disp([ mfilename ': Models_Edit: TODO' ])
  % create new Models after edition
  
function mifit_Models_Plot(varargin)
  % TODO
  disp([ mfilename ': Models_Plot: TODO' ])
  
function mifit_Models_Plot_Parameters(varargin)
  % TODO
  disp([ mfilename ': Models_Plot_Parameters: TODO uitable' ])

function mifit_Models_Export_Parameters(varargin)
  % * Export to file... (JSON, M, YAML, MAT...)
  disp([ mfilename ': Models_Export_Parameters: TODO' ])
  
function mifit_Models_Add_Expression(varargin)
  % * 4D TAS convolution        -> in Models Transformation/operations
  % * Powder average 4D -> 2D   -> in Models Transformation/operations
  % TODO: is this needed ?

% set optimizer configuration -> contextual dialogue in Model_Parameters uitable ?

% Optimizers menu **************************************************************

function mifit_Optimizers_Configure(varargin)
% Optimizers/Configure: open dialogue to change current optimizer parameters
  % change optimizer configuration parameters
  CurrentOptimizer = getappdata(mifit_fig,'CurrentOptimizer');
  mifit_disp([ '[Optimizers_Configure] Configure "' CurrentOptimizer '"' ]);
  options    = getappdata(mifit_fig, 'CurrentOptimizerConfig');
  config     = getappdata(mifit_fig, 'Preferences');
  o.FontSize = config.FontSize;
  o.Name     = [ mfilename ': Optimizer "' CurrentOptimizer '" configuration' ];
  NL = sprintf('\n');
  o.TooltipString = [ 'Most relevant "' CurrentOptimizer '" configuration items:' NL ...
    '* MaxFunEvals - Maximum number of function evaluations allowed [ positive integer ]' NL ...
    '* MaxIter - Maximum number of iterations allowed [ positive scalar ]' NL ...
    '* TolFun - Termination tolerance on the function value [ positive scalar ]' NL ...
    '* TolX - Termination tolerance on X [ positive scalar ]' ];
  options1 = structdlg(options, o);
  if isempty(options1), return; end
  % look for changes in new options...
  fields = fieldnames(options);
  for index=1:numel(fields)
    if ~isequal(options1.(fields{index}), options.(fields{index}))
      t1 = class2str(fields{index},options1.(fields{index})); t1 = strrep(t1, sprintf('\n'), '');
      t0 = class2str(fields{index},options.(fields{index}));  t0 = strrep(t0, sprintf('\n'), '');
      mifit_disp([ '[Optimizers_Configure] Assigned ' CurrentOptimizer ': ' ...
        t1 ' [was ' t0 ']' ]);
    end
  end
  options = options1;
  setappdata(mifit_fig,'CurrentOptimizerConfig', options);
  

function mifit_Optimizers_Criteria(varargin)
% Optimizers/Criteria: open dialogue to change the fit criteria
  
  crit = {'least_square = (|Signal-Model|/Error)^2 [non-robust]', ...
    'least_absolute = |Signal-Model|/Error [robust]',...
    'least_median = median(|Signal-Model|/Error) [robust]',...
    'least_mean = mean(|Signal-Model|/Error) [robust]',...
    'least_max = max(|Signal-Model|/Error) [non-robust]',...
    'least_rfactor = (|Signal-Model|/Error)^2/(Signal/Error)^2 [non-robust]',...
    'max_likelihood =  (|Signal-Model|/σ)^2 + 1/2 log(2πσ) [robust]', ...
    'max_corrcoef = 1-corrcoef(Signal, Model)'};
  selection = [];
  if isappdata(mifit_fig,'CurrentOptimizerCriteria')
    selection = find(strcmp(getappdata(mifit_fig,'CurrentOptimizerCriteria'), strtok(crit)),1);
  end
  if isempty(selection), selection = 1; end
  % show operator selection dialogue
  [selection, ok] = listdlg('ListString', crit, 'ListSize',[400 160], ...
    'Name',[ 'miFit: Select the Fit Criteria ti use' ], ...
      'PromptString', 'The Fit criteria is the quantity which is to be minimized. It is a function of the Data set Signal, and the Model value. It is then divided by the number of degrees of freedom.', ...
      'InitialValue',selection, 'SelectionMode','single');
  if isempty(selection) || ok ~= 1, return; end
  crit = strtok(crit{selection});
  setappdata(mifit_fig,'CurrentOptimizerCriteria', crit);

function selected = mifit_Optimizers_Set(varargin)
% [internal] mifit_Optimizers_Set: set current optimizer, and check its label in the menu
  % set the optimizer to use
  selected0 = getappdata(mifit_fig, 'CurrentOptimizer');
  selected  = selected0;
  if nargin && ~isempty(varargin{1})
    if ishandle(varargin{1})
      selected = get(varargin{1},'UserData');
    elseif ischar(varargin{1})
      selected = varargin{1};
    end
  end
  if isempty(selected), selected='fmin'; end  % default when not set
  if ~strcmp(selected, selected0) % when we change the optimizer, we reset its Configuration
    options = feval(selected, 'defaults');
    setappdata(mifit_fig,'CurrentOptimizerConfig', options);
  end
  mifit_disp([ '[Optimizer] Setting optimizer to "' selected '"' ]);
  setappdata(mifit_fig,'CurrentOptimizer', selected);

  % we get the CurrentOptimizer, and check it. uncheck the others
  hmodels = mifit_fig('Menu_Optimizers');
  hoptims = get(hmodels,'Children');
  names   = get(hoptims,'UserData');
  index   = strcmp(selected, names);
  set(hoptims,'Checked','off');
  set(hoptims(index),'Checked','on');
  
% Help menu *******************************************************************

function mifit_Help_Main(varargin)
% Help/Main: open web page for miFit
  % TODO
  disp([ mfilename ': Help_Main: TODO' ])
  
function mifit_Help_Loaders(varargin)
% Help/Loaders: open web page Loaders
  doc(iData,'Loaders');
  
function mifit_Help_Models(varargin)
% Help/Models: open web page Models
  doc(iData,'Models');
  
function mifit_Help_Optimizers(varargin)
% Help/Optimizers: open web page Optimizers
  doc(iData,'Optimizers');
  
function mifit_Help_Math(varargin)
% Help/Math: open web page Math
  doc(iData,'Math');
  
% List Data and Stack management ***********************************************

function mifit_List_Data_Files(varargin)
% [internal] mifit_List_Data_Files: called when clicking on the listbox
% TODO: support KeyPressFcn to e.g. 'p'=plot(log10(d)), 'f'=fit, 's'=save, 'o'=open, ...
% TODO: install a UIContextMenu in the Data sets list uicontrol ?
  if nargin && ishandle(varargin{1}), 
    obj=varargin{1}; 
  else obj=[]; end

  [d, index_selected] = mifit_List_Data_pull(); % get selected objects
  
  if isempty(index_selected) || numel(d) == 0
    setappdata(mifit_fig, 'CurrentDataSet', []);
    setappdata(mifit_fig, 'CurrentDataSetIndex', []);
  end

  if ~isempty(obj)
    if strcmp(get(gcbf,'SelectionType'), 'open')  % double-click
      mifit_Data_Plot(d);  % plot every time we change the data set selection
    end
  end
  
  % when a CurrentModel exists (previously selected from Models menu)
  % we assign that one to Data sets which have no Model defined
  if ~isempty(getappdata(mifit_fig, 'CurrentModel'))
    model = getappdata(mifit_fig, 'CurrentModel');
    D     = getappdata(mifit_fig, 'Data');  % all data sets
    model_assignments = 0;
    for index=index_selected(:)'
      if numel(D) > 1, this_d = D(index); else this_d = D; end
      if isempty(this_d) || ~isa(this_d,'iData'), continue; end
      previous_model = [];
      % get the Model stored in the Dataset (after fit)
      if isfield(this_d, 'Model')
        previous_model = get(this_d, 'Model');
      elseif ~isempty(findfield(this_d, 'Model'))
        previous_model = get(this_d, findfield(this_d, 'Model', 'cache first'));
      end
      if isempty(previous_model)  % store the CurrentModel in the Dataset
        this_d = setalias(this_d, 'Model', model);
        model_assignments = model_assignments+1;
      end
      if numel(D) > 1, D(index) = this_d; else D = this_d; end
    end
    if model_assignments
      setappdata(mifit_fig, 'Data', D);
      mifit_disp([ 'Assigning Model "' model.Name '" to ' num2str(model_assignments) ' Data set(s).' ]);
    end
  end
  
  % store the first selected Model so that its Model Parameters can be updated
  % in mifit_Models_View_Parameters
  if ~isempty(index_selected) && numel(d)
    d=d(1); index_selected=index_selected(1);
  end
  setappdata(mifit_fig, 'CurrentDataSet', d);
  setappdata(mifit_fig, 'CurrentDataSetIndex', index_selected);
  
  if ~isempty(mifit_fig('mifit_View_Parameters'))
    % trigger an update of the Parameter window when already opened
    mifit_Models_View_Parameters('update');
  end

% ******************************************************************************


  

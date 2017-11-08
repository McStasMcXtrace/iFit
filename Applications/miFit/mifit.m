function varargout = mifit(varargin)
% miFit: a user interface to iFit
%
% data  = mifit('data')      retrieves all data sets from the interface
% models= mifit('models')    retrieves user models from the interface
% data  = mifit('pull')      retrieves only selected data sets
% mifit('push', datasets)   replace existing data sets and append new ones
% mifit('push', models)     replace existing models and append new ones
% mifit('filename')         imports the file into a new Data set/Model
% mifit(iData_object)       add the iData object into the interface Data stack
% mifit(iFunc_object)       add the iFunc Model into the interface Models menu
% config = mifit('config')  retrieves the miFit configuration
% mifit('exit')             closes miFit
% mifit('reset')            reset miFit to its Factory settings
% mifit('handles')          retrieves the miFit GUI handles
%
% Version: $Date$
% (c) E.Farhi, ILL. License: EUPL.

%  function config = mifit_Preferences_Load
%  function mifit_Preferences_Save(config)
%  function mifit_File_New(handle)
%  function mifit_File_Open(handle)
%  function mifit_File_Save(varargin)
%  function mifit_File_Saveas(varargin)
%  function mifit_File_Print(varargin)
%  function mifit_File_Preferences(varargin)
%  function mifit_Preferences_Apply
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
    if isempty(fig) && nargin == 1
      % interface does not exist yet. Handle high priority commands.
      if strcmp(varargin{1},'File_Exit') || strcmpi(varargin{1},'exit'), return;
      elseif strcmp(varargin{1},'File_Reset_Factory') ...
      || strcmp(varargin{1},'reset') || strcmp(varargin{1},'factory')
        mifit_disp([ '[RESET] Factory defaults. ' ]);
        file = fullfile(prefdir, [ mfilename '.ini' ]);
        if ~isempty(dir(file)), delete(file); end
        file = fullfile(prefdir, [ mfilename '.mat' ]);
        if ~isempty(dir(file)), delete(file); end
        return
      end
    end
    
    if isempty(fig) || ~ishandle(fig)
      fig = feval([ mfilename '_OpeningFcn' ]);
    else set(fig, 'NextPlot','new','HandleVisibility','callback');
    end
    out = fig;
    if ~isempty(varargin)
        if nargin == 2 && isa(varargin{1}, 'dndcontrol') && isstruct(varargin{2})
          % ** Drag-n-Drop: called as a callback/event from Drag-n-Drop (import external file/text)
          evt = varargin{2};
          switch evt.DropType
          case 'file'
              for n = 1:numel(evt.Data)
                  mifit(deblank(evt.Data{n}));
              end
          case 'string'
              lines = textscan(evt.Data, '%s','Delimiter',sprintf('\n'));
              mifit(deblank(lines{1}));
          end
        elseif ischar(varargin{1}) && isempty(dir(varargin{1})) ...
          && ~any(strncmp(varargin{1}, {'http:','ftp:/','file:','https'}, 5)) % a function/action to call ?
          % ** ACTION: callback with varargin{1} == 'action'
          action = varargin{1};
          if strcmpi(action,'identify'), varargout{1} = []; return; end
          if any(strcmpi(action,{'data'}))
            % get the full data list
            out = getappdata(fig, 'Data');
          elseif any(strcmpi(action,{'models'}))
            out = mifit_Models_GetList();
          elseif any(strcmpi(action,{'pull','selection'}))
            % get the selected data sets
            out = mifit_List_Data_pull();
          elseif any(strcmpi(action,{'config'}))
            % get the configuration
            out = getappdata(mifit_fig, 'Preferences');
          elseif any(strcmpi(action,{'handles'}))
            % get the configuration
            out = mifit_fig('handles');
          elseif any(strcmpi(action,{'push','replace','merge'})) && nargin >= 2
            mifit_List_Data_push(varargin(2:end), 'replace');
            return
          elseif any(strcmpi(action,{'exit','quit'}))
            mifit_File_Exit;
            return
          elseif any(strcmpi(action,{'reset','factory'}))
            mifit_File_Reset_Factory(action);
            return
          elseif ~isempty(action)
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
                  config1 = mifit_Preferences_Load(file);
                catch
                  config1 = [];
                end
                if ~isempty(config1)
                  mifit_Preferences_Apply(config1);
                  % and merge with existing config
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

function mifit_File_Eval(handle)
  mifit_disp([ '[File_Eval] Type Matlab code in the TextEdit window, and ' ...
      'evaluate it from the File menu. You can use the normal Matlab syntax, ' ...
      'including tests (if ... else ... end), loops (for ... end), but ' ...
      'commands should all end with ";". To display some result, use disp(variable).']);
  TextEdit;
  
function mifit_File_Print(varargin)
% File/Print: print the interface. 
% generate an HTML report and display in browser for printing.
  d=mifit_List_Data_pull(); % get selected objects
  d = { d{:} };
  if all(isempty(d)), return; end
  dirname  = tempname;
  filename = fullfile(dirname,'index.html');
  mkdir(dirname);
  mifit_disp([ '[File_Print] Exporting Data sets to HTML ' filename ' for printing...' ]);
  save(d, filename, 'html data');
  webbrowser(filename,'system');  % tries to open with the system browser

  
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
  [d, index_selected]=mifit_List_Data_pull(); % get selected objects
  if ~isempty(index_selected)
    d=d{1};
    setappdata(mifit_fig, 'CurrentDataSet', d);
    setappdata(mifit_fig, 'CurrentDataSetIndex', index_selected);
    if ~isempty(mifit_fig('mifit_View_Parameters'))
      % trigger an update of the Parameter window when already opened
      mifit_Models_View_Parameters('update');
    end
  end

function mifit_Edit_Cut(varargin)
% Edit/Cut: get the selected indices in the List, copy these elements to the clipboard
% and delete the elements. Update the History (in Delete).
  mifit_Edit_Copy(varargin{:});
  mifit_Edit_Delete(varargin{:});



function mifit_Edit_Duplicate(varargin)
% Edit/Duplicate: copy and paste selected. Update the history (in Paste).
  d=mifit_List_Data_pull();
  for index=1:numel(d)
    d{index} = copyobj(d{index});
  end
  mifit_List_Data_push(d);

  
% Data menu ********************************************************************

  
function mifit_Data_PlotAs_Overlay(varargin)
% Data/Plot: Overlay
  mifit_Data_PlotAs('', varargin{:});
  
function mifit_Data_PlotAs_Contour(varargin)
% Data/Plot: Contour
  mifit_Data_PlotAs('contour3 light transparent grid tight replace', varargin{:});
  
function mifit_Data_PlotAs_Waterfall(varargin)
% Data/Plot: Waterfall
  mifit_Data_PlotAs('waterfall light transparent grid tight replace', varargin{:});
  
function mifit_Data_PlotAs_Surface(varargin)
% Data/Plot: Surface
  mifit_Data_PlotAs('surfc light transparent grid tight replace', varargin{:});
  
function mifit_Data_PlotAs_ScatterPlot(varargin)
% Data/Plot: Scatter
  mifit_Data_PlotAs('scatter3 light transparent grid tight replace', varargin{:});
  
function mifit_Data_PlotAs_Plot3(varargin)
% Data/Plot: Plot3
  mifit_Data_PlotAs('plot3 grid tight replace', varargin{:});
  
function mifit_Data_PlotAs_Slice(varargin)
% Data/Plot: Slice (2D/3D)
  mifit_Data_PlotAs('slice', varargin{:});

function mifit_Data_Saveas(varargin)
% Data/Saveas: export selected data sets
  d = mifit_List_Data_pull;
  if numel(d)
    save([ d{:} ], 'gui');
  end
  
function mifit_Data_Table(varargin)
% Data/view Table: view data sets as tables
  d = mifit_List_Data_pull;
  config = getappdata(mifit_fig, 'Preferences');
  
  for index=1:numel(d)
    handle = edit(d{index}, 'editable');
    set(handle, 'FontSize', config.FontSize);
    set(handle, 'DeleteFcn', @mifit);
  end

function mifit_Data_View(varargin)
% Data/view Source: edit source file when exists
  d = mifit_List_Data_pull;
  for index=1:numel(d)  
    if ~isempty(d{index}.Source) && ~isdir(d{index}.Source)
      try
        edit(d{index}.Source)
      end
    end
  end
  
function mifit_Data_Properties(varargin)
% Data/Properties: display dialogue about data sets
  d = mifit_List_Data_pull;
  help([ d{:} ]);  % this is very basic. No edit there.


  
function mifit_Data_History(varargin)
% Data/History: show Data sets command history
  d = mifit_List_Data_pull();
  for index=1:numel(d)
    [c,fig]=commandhistory(d{index});
  end
  
% Models menu ******************************************************************



function mifit_Models_Load(varargin)
% Models/Import: load new Models
% we pop-up a file selector, read as iFunc, and store it in miFit.
  models = load(iFunc, '');
  if isempty(models), return; end
  mifit_disp([ '[Models] Loading ' num2str(numel(models)) ' new models:' ]);
  mifit_disp(display(models));
  mifit(models);

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


% Help menu *******************************************************************

function mifit_Help_Main(varargin)
% Help/Main: open web page for miFit
  % TODO
  doc(iData,'miFit');
  
function mifit_Help_Loaders(varargin)
% Help/Loaders: open web page Loaders
  disp(' ');
  data1 = iLoad('formats');
  disp(' ');
  data2 = save(iData, 'formats');
  disp(' ');
  data3 = save(iFunc, 'formats');
  disp(' ');
  if isfield(data1, 'loaders')
    mifit_disp([ '[Help] Data loaders:   ' num2str(numel(data1.loaders)) ' formats' ]);
  end
  
  if size(data2,1) > 1
    mifit_disp([ '[Help] Data exporters: ' num2str(size(data2,1)) ' formats' ]);
  end
  
  if size(data3,1) > 1
    mifit_disp([ '[Help] Model exporters:' num2str(size(data3,1)) ' formats' ]);
  end
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

% ******************************************************************************


  

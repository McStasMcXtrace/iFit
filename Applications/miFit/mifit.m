function out = mifit(varargin)


% notes:
% saveas: display dialogue to select what to export:
% * data sets
% * parameters
% * fit model

% plot: plot model/dataset/parameters ? -> preferences

% set optimizer configuration -> contextual dialogue

% treat data set with MView type operations

% re-asign data set signal, axes, ... to aliases/new ones

% will add toolbar instead of uicontrols for fast processing

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
    end

    if ~isempty(varargin) 
        if ischar(varargin{1})          % a function/action to call ?
            % callback with varargin{1} == 'action'
            if strcmpi(varargin{1},'identify'), return; end
            try
              feval(varargin{:});
            catch ME
              disp([ mfilename ': Unknown action "' varargin{1} '"' ])
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
                List_Data_push(d);
            end
        end
    end
    
% -------------------------------------------------------------------------
% Main window ID and creation
% -------------------------------------------------------------------------

function f=mifit_fig(tag)
% search for a given Tag in Application or main Figure if ommitted.
  persistent fig handles
  
  if ~ishandle(fig), fig=[]; end
  if isempty(fig)
    fig = findall(0, 'Tag','MiFit');
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
    disp([ mfilename ': Welcome to MiFit !' ])
    disp(datestr(now))
    fig = openfig(mfilename);
    
    % test if any of the Models/Optimizers menu is empty
    % if not, return (no need to build them)
    hmodels = mifit_fig('Menu_Models');
    hoptim  = mifit_fig('Menu_Optimizers');
    if isempty(hmodels) || isempty(hoptim) ...
            || ~isempty([ get(hmodels,'Children') ; get(hoptim,'Children') ]), return; end
    
    % Display welcome dialog during menu build
    h = msgbox('Welcome to MiFit !','MiFit: Starting','help');
    
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
    
    % create the AppData
    % Data Stack
    setappdata(fig, 'Data', []);
    
    % close welcome image
    delete(h);
end

% -------------------------------------------------------------------------
% Callbacks
% -------------------------------------------------------------------------

function File_New(handle)
% File/New menu item
  d = iData(zeros(5)); % create an empty Data set;
  disp([ mfilename ': Editing an empty Data set. Close the window to retrieve its content' ]);
  disp( '    as a new Data set into MiFit. Use the Contextual menu.');
  handle = edit(d, 'editable');
  % set a DeletedFcn so that the content can be retrieved into MiFit when
  % closing.
  set(handle, 'DeleteFcn', @mifit);
  
function File_Open(handle)
  d = iData('');
  % push that data onto the List
  
function List_Data_push(d)
% put a new data set at the end of the stack
  fig = mifit_fig;
  Data = getappdata(fig, 'Data');
  
  % update AppData Stack
  if numel(d) > 1, d = d(:); end
  Data = [ Data ; d ];  % a columns of iData set
  setappdata(fig, 'Data', d);
  
  % update the List labels by appending the Name at the end
  hObject        = mifit_fig('List_Data_Files');
  list           = get(hObject,'String');
  index_selected = get(hObject,'Value');
  for index=1:numel(d)
      list{end+1} = char(d);
      index_selected(end+1) = 1;
  end
  set(hObject,'String', list, 'Value', index_selected);
  

function d=List_Data_Files(hObject)
% get the selected Data List
  if nargin == 0 || isempty(hObject)
      hObject = mifit_fig('List_Data_Files');
  end
  d = [];

  index_selected = get(hObject,'Value');
  list = get(hObject,'String');
  if isempty(list), return; end
  
  item_selected = list{index_selected};
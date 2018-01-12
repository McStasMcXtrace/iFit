function fig = mifit_OpeningFcn
% [internal] mifit_OpeningFcn: This function creates the main mifit window and returns its ID

fig = mifit_fig();
if isempty(fig) || ~ishandle(fig)
    % create the main figure
    mifit_disp('[Init] Welcome to miFit ! ******************************************************')
    fig = openfig('mifit');
    set(fig, 'HandleVisibility','callback','NumberTitle','off','NextPlot','new');
    
    % get Models/Optimizers menu handles
    hmodels = mifit_fig('Menu_Model');
    hoptim  = mifit_fig('Menu_Optimizers');
    
    % load Preferences
    mifit_Preferences_Load;
    config=mifit_Preferences_Apply;

    % create the AppData default values
    setappdata(fig, 'Data',    {});
    setappdata(fig, 'History', {});
    setappdata(fig, 'Models',  {});
    setappdata(fig, 'CurrentOptimizer',   []);
    setappdata(fig, 'CurrentOptimizerCriteria',[]);
    setappdata(fig, 'CurrentOptimizerConfig',[]);
    setappdata(fig, 'CurrentModel',       []);
    setappdata(fig, 'CurrentModelHandle', []);
    setappdata(fig, 'CurrentDataSet',     []);
    setappdata(fig, 'CurrentDataSetIndex',[]);
    setappdata(fig, 'CurrentDataSetHandle',[]);
    
    set(fig,'Pointer','watch');
    
    % Display welcome dialog during menu build
    h = mifit_Help_About(fig);
    mifit_disp([ version(iData) sprintf(' Visit <http://ifit.mccode.org>') ]);

    % get the list of Models and Optimizers
    models = []; d=[];
    file = fullfile(prefdir, [ 'mifit' '.mat' ]);
    if ~isempty(dir(file))
      try
        d = load(file);
        if isfield(d, 'Models')
          fd = dir(file);
          mifit_disp([ '[Init] Loading Model library from ' file ' [size: ' num2str(round(fd.bytes/1024)) ' kb]' ]);
          models = d.Models;  % contains the callback and the label members
        end
      end
    end
    if isempty(models)
      mifit_disp([ '[Init] Building Optimizer and Model library from ' file '. Be patient (only once)...' ]);
      [optimizers,models,filenames] = fits(iFunc);
      % build Models.callback and Models.label. We will instantiate at callback.
      modelsF = cell(1, numel(models));
      for index=1:numel(models)
        this = struct('callback',filenames{index}, ...
          'label',get(models(index), 'Name'), 'object',models(index));
        modelsF{index}=this;
      end
      models = modelsF;
    else
      optimizers = fits(iFunc);
    end
    
    % fill Models menu
    mifit_Models_Add_Entry(models);
    
    % fill Optimizers menu
    if ~isempty(optimizers) && iscell(optimizers)
        mifit_disp([ '[Init] Initializing ' num2str(numel(optimizers)) ' Optimizers ...' ]);
        separator = 'on';
        for f=optimizers
            % each optimizer is given with its function name. We request
            % 'defaults' and display its name
            
            o=feval(f{1},'defaults');
            if isfield(o, 'algorithm') && ~isempty(o.algorithm)
                algorithm = o.algorithm;
            else
                algorithm = f{1};
            end
            algorithm = strtrim(algorithm);
            if ~isempty(algorithm)
              % search if this entry already exists
              item_handle = findobj(hoptim, 'Type','uimenu', 'Label', algorithm);
              if isempty(item_handle)
                uimenu(hoptim, 'Label', algorithm, 'UserData', f{1}, ...
                  'Callback', 'mifit(''Optimizers_Set'',gcbo);', ...
                  'separator', separator);
                separator = 'off';
              end
            end
        end
    end
    
    % assign the saved CurrentOptimizer and other saved stuff
    if ~isfield(d, 'CurrentOptimizer'), d.CurrentOptimizer = []; end
    mifit('Optimizers_Set',d.CurrentOptimizer);
    setappdata(fig, 'Optimizers',optimizers);
    if isfield(d,'CurrentOptimizerConfig'),   setappdata(fig, 'CurrentOptimizerConfig',   d.CurrentOptimizerConfig); end
    if isfield(d,'CurrentOptimizerCriteria'), setappdata(fig, 'CurrentOptimizerCriteria', d.CurrentOptimizerCriteria); end
    if isfield(d,'CurrentModel'),             setappdata(fig, 'CurrentModel',             d.CurrentModel); end
    
    % add predefined configurations to the File menu
    % search for item Configurations
    handle = findobj(fig, 'Tag','File_Configurations');
    path_config = fullfile(fileparts(which(mfilename)),'..','configurations');
    % insert Configurations items
    for conf = dir(path_config)'
      if ~conf.isdir && conf.name(end) ~= '~' % not dir nor temp file
        uimenu(handle, 'Label', conf.name, ...
                'Callback', [ 'mifit(''' fullfile(path_config, conf.name) ''');' ]);
      end
    end
    
    % create the AppData Data Stack
    % Load the previous Data sets containing Model Parameters (when a fit was performed)
    if ~isempty(d)
      if isfield(d, 'Data')
        mifit_disp([ '[Init] Loading Data sets from ' file ]);
        mifit_List_Data_push(d.Data);
      end
    end
    
    % resize main panel
    handle = findobj(fig, 'Tag','Panel_DataSets');
    set(handle, 'Units','normalized', 'Position',[0.002 0 0.99 0.99]);
    
    % start Logging
    file = fullfile(prefdir, [ 'mifit' '.log' ]);
    mifit_disp([ '[Init] Log file is ' file ]);

    % initialize Java hooks and activate Drag-n-Drop from external source (files, text)
    if ~(exist('org.yaml.snakeyaml.Yaml','class') == 8)
      javaaddpath(YAML.JARFILE);
    end
    hObject = mifit_fig('List_Data_Files');
    if ~(exist('MLDropTarget','class') == 8)
      dndcontrol.initJava;
    end
    dnd = dndcontrol(hObject,@mifit,@mifit);
    setappdata(fig, 'DnDControl', dnd);

    % activate Drag-n-Drop from the List to other Matlab windows
    % set(fig,'windowbuttonupfcn',  'disp(''up in:''); get(0,''pointerwindow'')')
    % set(fig,'windowbuttondownfcn','disp(''down in:''); get(0,''pointerwindow'')')

    % close welcome image
    if ishandle(h)
      delete(h);
    end
    
    set(fig,'Pointer','arrow');

end

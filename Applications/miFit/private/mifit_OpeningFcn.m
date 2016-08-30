function fig = mifit_OpeningFcn
% This function creates the main window and returns its ID

fig = mifit_fig();
if isempty(fig) || ~ishandle(fig)
    % create the main figure
    mifit_disp('[Init] Welcome to miFit ! ********************************************')
    fig = openfig('mifit');
    
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
    file = fullfile(prefdir, [ 'mifit' '.mat' ]);
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
      modelsF = cell(1, numel(models));
      for index=1:numel(models)
        this = struct('callback',filenames{index},'label',get(models(index), 'Name'),'object',models(index));
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
    setappdata(fig, 'Optimizers',optimizers);
    
    % Load the previous Data sets containing Model Parameters (when a fit was performed)
    file = fullfile(prefdir, [ 'mifit' '.mat' ]);
    if ~isempty(dir(file))
      try
        d = load(file);
        if isfield(d, 'Data')
          mifit_disp([ '[Init] Loading Data sets from ' file ]);
          mifit_List_Data_push(d.Data);
        end
      end
    end
    file = fullfile(prefdir, [ 'mifit' '.log' ]);
    mifit_disp([ '[Init] Log file is ' file ]);

    % close welcome image
    delete(h);

end

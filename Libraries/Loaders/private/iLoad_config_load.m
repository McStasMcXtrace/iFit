function config = iLoad_config_load
% iLoad_config_load: load the configuration and format customization

  loaders      = {};
  config       = [];
  % read user list of loaders which is a cell of format descriptions
  if exist(fullfile(prefdir, 'iLoad.ini'), 'file')
    % there is an iLoad_ini in the Matlab preferences directory: read it
    configfile = fullfile(prefdir, 'iLoad.ini');
    content    = fileread(configfile);
    % evaluate content of file
    try
      eval(content(:)'); % this makes a 'config' variable
    end
    if iscell(config)
      loaders = config; config=[];
    elseif isstruct(config) && isfield(config, 'loaders')
      loaders = config.loaders;
    end
    if ~isempty(loaders)
      config = [];  % reset struct
      config.FileName= configfile;
      disp([ '% Loaded iLoad format descriptions from ' config.FileName ]);
    end
  end
  % append default iLoad_ini configuration
  if exist('iLoad_ini', 'file')
    defaults = iLoad_ini;
    if isstruct(defaults) && isfield(defaults, 'loaders')
      loaders = [ loaders(:) ; defaults.loaders(:) ];
    elseif iscell(defaults)
      loaders = [ loaders(:) ; defaults(:) ];
    else defaults = [];
    end
    if ~isfield(config, 'FileName'), config.FileName= which('iLoad_ini'); end
  end
  
  % read the output of Loaders/read_* functions ('identify', return struct)
  p = fileparts(which('iLoad'));
  p = dir(fullfile(p, 'read_*'));
  for index=1:numel(p)
    [~,this] = fileparts(p(index).name);  % skip M extension
    try
      id = feval(this, 'identify');
    catch ME
      disp([ mfilename ': WARNING: ' this '(''identify''): invalid loader information.' ]);
      id = [];
    end
    % the output can be an array: append elements to loaders list
    for jj=1:numel(id)
      if iscell(id)
        this = id{jj};
      elseif isstruct(id)
        this = id(jj);
      else this=[];
      end
      if isstruct(this) && isfield(this, 'name') && isfield(this, 'method')
        loaders{end+1} = this;
      end
    end
  end
  
  % ADD default loaders: method, ext, name, options
  % default importers, when no user specification is given. 
  % format = { method, extension, name, {options, patterns, postprocess} }
  formats = {...
    { 'csvread', 'csv', 'Comma Separated Values (.csv)','',''}, ...
    { 'dlmread', 'dlm', 'Numerical single block','',''}, ...
    { 'xlsread', 'xls', 'Microsoft Excel (first spreadsheet, .xls)','',nan}, ...
    { 'load',    'mat', 'Matlab workspace (.mat)',nan,'','openhdf'}, ...
    { 'importdata','',  'Matlab importer','',nan}, ...
  };
  for index=1:length(formats) % the default loaders are addded after the INI file
    format = cell(1,6);
    [format(1:numel(formats{index}))] = deal(formats{index});

    loader.method     = format{1};
    loader.extension  = format{2};
    loader.name       = format{3};
    loader.options    = format{4};
    loader.patterns   = format{5};  % patterns=nan for binary formats; '' for text without patrern
    loader.postprocess= format{6};
    loaders{end+1} = loader;
  end
  
  % sort loaders. Only get the first entry for each (remove duplicates).
  % the loaders sort must preserve the initial order. Collect patterns for sorting.
  names = {}; patterns={}; loaders=loaders(:);
  for index=1:numel(loaders)
    loader = loaders{index};
    % remove invalid loaders, only add when not already registered in 'names'
    if isfield(loader, 'method') && isfield(loader, 'name') && ...
      (isa(loader.method, 'function_handle') || exist(loader.method,'file')) && ...
      ~any(strcmp(loader.name, names))
        names{end+1}    = loader.name;
        patterns{end+1} = nan; % default is no pattern, assumed binary
        jj = numel(names);
        if isfield(loader, 'patterns')
          if ischar(loader.patterns)
            patterns{jj} = { loader.patterns };
          elseif iscellstr(loader.patterns)
            patterns{jj} = loader.patterns;
          end
        end
    else names{end+1} = ''; patterns{end+1} = nan; end
  end
  index   = find(~cellfun(@isempty, names));
  loaders = loaders(index); % first valid and unique
  patterns= patterns(index);
  names   = names(index);
 
  % then put first 'text' formats with patterns (longest first), then 'text' without patterns
  % then 'binary' ones (patterns=nan)
  elem    = cellfun(@numel, patterns) + cellfun(@iscellstr, patterns);
  index1  = find(elem > 1);
  loaders1 = loaders(index1); % get formats with patterns
  patterns1= patterns(index1);
  loaders(index1) = [];       % remove these form formats, keep others in original order
  % sort loaders with patterns, decreasing specificity
  [elem1,index1]   = sort(cellfun(@numel, patterns1),2,'descend');
  loaders1 = loaders1(index1);
  loaders  = [ loaders1 ; loaders ];  % patterns up-front
  
  % check if other configuration fields are present, else defaults
  if ~isfield(config, 'UseSystemDialogs'), config.UseSystemDialogs = 'no'; end
  if ~isfield(config, 'FileName'),         config.FileName = ''; end
  if ~isfield(config, 'MeX'), config.MeX = []; end
  if ~ischar(config.MeX) && ~isempty(config.MeX)
    if ~isfield(config.MeX, 'looktxt')
      if ispc || ismac, config.MeX.looktxt = 'yes'; 
      else              config.MeX.looktxt = 'no'; end % Linux: avoid MeX which may be unstable (SEGV)
    end
    if ~isfield(config.MeX, 'cif2hkl'),      config.MeX.cif2hkl = 'yes'; end
  end
  
  config.loaders = loaders; % updated list of loaders
  
end % iLoad_config_load

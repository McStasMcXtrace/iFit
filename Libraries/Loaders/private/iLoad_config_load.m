function config = iLoad_config_load
% iLoad_config_load: load the configuration and format customization

  persistent cached_config

  if ~isempty(cached_config)
    config = cached_config;
    return
  end
  
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
    { 'csvread', 'csv', 'Comma Separated Values (.csv)',            '',''}, ...
    { 'dlmread', 'dlm', 'Numerical single block',                   '',''}, ...
    { 'xlsread', 'xls', 'Microsoft Excel (first spreadsheet, .xls)','',[]}, ...
    { 'load',    'mat', 'Matlab workspace (.mat)',                  '',[],'openhdf'}, ...
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
  
  % make sure loaders contain the proper fields
  s = {};
  for index=1:numel(loaders)
    this = iLoad_check_loader(loaders{index});
    if ~isempty(this), s{end+1} = this; end
  end
  loaders = s;
  
  % get unique entries (this also sorts in alpha)
  [~, index] = unique(cellfun(@(c)getfield(c, 'name'), loaders,'UniformOutput',false));
  loaders = loaders(index);
  
  % add the default Matlab importer, but must be last
  loaders{end+1} = iLoad_check_loader( ...
    struct('name','Matlab importer','method','importdata'));
  
  % we sort loaders with highest pattern counts and text on top
  patterns   = cellfun(@(c)getfield(c, 'patterns'), loaders,'UniformOutput',false);
  istext     = cellfun(@(c)getfield(c, 'istext'), loaders);
  ext        = cellfun(@(c)getfield(c, 'extension'), loaders,'UniformOutput',false);
  count      = cellfun(@(c)numel(char(c)), patterns)+istext+cellfun(@numel, patterns)+(~cellfun(@isempty,ext));
  [~,index]  = sort(count, 2, 'descend');
  loaders    = loaders(index);
  
  % check if other configuration fields are present, else defaults
  if ~isfield(config, 'verbosity'), config.verbosity = 1; end
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

  cached_config = config;
  
end % iLoad_config_load

function loader=iLoad_check_loader(loader)
% check that a given loader (struct) is standadised
  
  if ~isstruct(loader) || ~isfield(loader, 'method') || ~isfield(loader, 'name')
    loader=[]; return; 
  end
  if isempty(loader.method) || isempty(loader.name)
    loader=[]; return; 
  end
  if ~isa(loader.method, 'function_handle') && ~exist(loader.method,'file')
    loader=[]; return; 
  end
  
  if ~isfield(loader,'patterns'),  loader.patterns =[]; end
  if ~isfield(loader,'extension'), loader.extension=[]; end
  if ~isfield(loader,'options'),   loader.options  =[]; end
  if ~isfield(loader,'postprocess'),loader.postprocess  =[]; end
  
  if ischar(loader.extension)
    loader.extension = cellstr(loader.extension);
  end
  
  % check if 'text'
  loader.istext = false;
  if ischar(loader.patterns) || (~isempty(loader.patterns) && iscellstr(loader.patterns))
    loader.patterns = cellstr(loader.patterns);
    loader.istext = true;
  end
  if ~isempty(strfind(lower(loader.name),'text'))
    loader.istext = true;
  end
end

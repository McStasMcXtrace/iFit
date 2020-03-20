function [data, loader] = iLoad_import(filename, loader, config, varargin)
% function to import single data with given method(s)

  data = []; isbinary=0;
  if isfield(config, 'verbosity')
    verbose = config.verbosity;
  else
    verbose = 0; % set this to 1 to get more output for debugging
  end
  
  if isempty(dir(filename))
    loader = 'Failed to load file (does not exist)';
    return
  end
  
  % check loader
  if isempty(loader), loader='auto'; end
  if strcmp(loader, 'auto')
    [loader, isbinary] = iLoad_loader_auto(filename, config);
  elseif strcmp(loader, 'gui')
    [dummy, filename_short, ext] = fileparts(filename);
    filename_short = [ filename_short ext];
    [loader, isbinary] = iLoad_loader_auto(filename, config);
    loader_names=[loader{:}];
    tmp         = cell(size(loader_names)); tmp(:)={' - '};
    loader_names= strcat({loader_names.name}, tmp,{loader_names.method});
    [loader_index,ok]=listdlg(...
      'PromptString',...
        {'Select suitable import methods',['to import file ' filename_short ]}, ...
      'SelectionMode','Multiple',...
      'ListString', loader_names, ...
      'ListSize', [300 160], ...
      'Name', ['Loader for ' filename_short ]);
    if isempty(loader_index), loader=[]; return; end
    loader=loader(loader_index);
  elseif ischar(loader) || iscellstr(loader)
    % test if loader is the user name of a function
    loader = iLoad_config_find(loader);
  end
  
  % handle multiple loaders (cell or struct array)
  if (iscell(loader) || isstruct(loader)) && length(loader) > 1
    loader=loader(:);
    recompile = 0;
    for index=1:length(loader)
      if iscell(loader), this_loader = loader{index};
      else this_loader = loader(index); end
      try
        data = iLoad_import(filename, this_loader, config, varargin{:});
      catch ME
        l=ME;
        warning(l.message);
        warning(getReport(ME))
        [dummy, name_short, ext] = fileparts(filename);
        if verbose
          fprintf(1, 'iLoad: Failed to import file %s with method %s (%s). Ignoring.\n', name_short, this_loader.name, char(this_loader.method));
        end
        data = [];
        if strcmp(l.identifier, 'MATLAB:nomem') || any(strncmpi(l.message, 'out of memory',length('out of memory')))
          fprintf(1,'iLoad: Not enough memory. Skipping import of this file.\n');
          break;
        end
      end
      if ~isempty(data)
        loader = this_loader;
        return;
      end
    end % for
    loader = 'Failed to load file (all known methods failed)';
    return; % all methods tried, none effective
  end % if iscell
  if iscell(loader) & length(loader) == 1
    loader = loader{1};
  end

  % handle single char loaders (IMPORT takes place HERE) =======================
  if ischar(loader)
    tmp=loader; clear loader;
    loader.method = tmp; loader.options='';
  end
  if ~isfield(loader,'method'), return; end
  if ~isfield(loader,'name'), loader.name = loader.method; end
  if ~isfield(loader,'extension'), loader.extension = ''; end
  if ~isfield(loader,'options'), loader.options = ''; end
  if ~isfield(loader,'postprocess'), loader.postprocess = ''; end
  if isempty(loader.method), return; end
  % skip SVN/GIT/CVS files
  f = lower(filename);
  if any([ strfind(f, [ '.svn' filesep ]) ...
           strfind(f, [ '.git' filesep ]) ...
           strfind(f, [ '.cvs' filesep ]) ])
    loader = [];
    return
  end

  if verbose
    fprintf(1, 'iLoad: Importing file %s with method %s (%s)\n', filename, loader.name, loader.method);
  end
  % we select the calling syntax which matches the number of input arguments
  if iscellstr(loader.options)
    varg = { filename, loader.options{:}, varargin{:} };
  elseif ischar(loader.options) && ~isempty(loader.options) && nargin(loader.method) == 1
    if ~isempty(loader.options), loader.options = [ ' ' loader.options ]; end
    varg = { [ filename loader.options ], varargin{:} };
  elseif ischar(loader.options) && ~isempty(loader.options)
    varg = { filename, loader.options, varargin{:} };
  else
    varg = { filename, varargin{:} };
  end
  % reduce the number of input arguments to the one expected
  if ~any(strcmp(loader.method, {'text','read_anytext','looktxt'})) && nargin(loader.method) > 0
    narg = min(length(varg), nargin(loader.method));
    varg = varg(1:narg);
  end
  % avoid calling text importer with binary
  if any(strcmp(loader.method, {'text','read_anytext','looktxt'})) && isbinary, return; end
  
  % call the loader
  data = feval(loader.method, varg{:});
  if isempty(data), return; end

  % special test to avoid reading binary file with read_anytext
  data = iLoad_loader_check(filename, data, loader);
  
end % iLoad_import

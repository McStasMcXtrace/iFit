function [data, format] = iLoad(filename, loader, varargin)
% [data, loader] = iLoad(file, loader, ...) imports any data into Matlab. 
%
% DESCRIPTION:
% ------------
% Imports any data into Matlab, returned as a structure or cell.
%
%   Default supported formats include: any text based including CSV, Lotus1-2-3, SUN sound, 
%     WAV sound, AVI movie, NetCDF, FITS, XLS, BMP GIF JPEG TIFF PNG ICO images,
%     HDF4, HDF5, MAT workspace, XML, CDF, JSON, YAML, IDL, Python Numpy
%   Other specialized formats include: McStas, ILL, SPEC, ISIS/SPE, INX, EDF, Mantid.
%     SIF, MCCD/TIFF, ADSC, CBF, Analyze, NifTI, STL,PLY,OFF, CIF/CFL,
%     EZD/CCP4, Bruker Varian and JEOL NMR formats, Bruker OPUS, LabView LVM and TDMS,
%     Agilent and Thermo Finnigan MS, Quantum Design VMS, ENDF/ACE, Agilent SDF
%   Compressed files are also supported, with on-the-fly extraction (zip, gz, tar, Z).
%
%   To get a full list of supported formats, use at prompt: iLoad formats
%
%   Distant files are supported through e.g. URLs such as 
%     file://, ftp:// http:// and https://
%   File names may end with an internal anchor reference '#token", as usual in HTML 
%     links, in which case the members matching the token are returned.
%
% Typical usage:
%   iLoad('file'); 
%   iLoad('file.zip')
%   iLoad('file#token');
%   iLoad('http://path/name'); 
%
% Additional remarks:
% --------------------
% The file formats cache and MeX files can be rebuilt/checked with
%   iLoad force         update loaders cache and check importers
%   iLoad compile       force compilation of external importers (MeX)
%   iLoad save          save the cache to the user Preferences
%   iLoad config        get the loaders list as a configuration structure
%   iLoad formats       list of all supported file formats
%   iLoad FMT formats   list of file formats matching FMT
%
% The iLoad_ini configuration file can be loaded and saved in the Preference 
% directory using 
%   [config, configfile] = iLoad('load config').
%   [config, configfile] = iLoad(config,'save config')
%
% ------------------------------------------------------------------------------
% syntax:
%   [data, loader] = iLoad(file, 'loader', ...)
%
% input arguments:
%   file:   file name, or cell of file names, or any Matlab variable, or a URL
%             or an empty string (then pops'up a file selector)
%   loader: a function name to use as import routine, OR a structure with:
%             loader = 'auto' (default) 
%             loader = 'gui'  (asks for the format to use)
%             loader.method = 'function name'
%             loader.options= string of options to catenate after file name
%                          OR cell of options to use as additional arguments
%                             to the method
%           leave it empty ('') for automatic setting.
%   additional arguments are passed to the import routine.
%
% output variables:
%   data:   a single structure containing file data, or a cell of structures
%   loader: the loader that was used for importation, or a cell of loaders.
%
% Example: a=iLoad(fullfile(ifitpath, 'Data')); iscell(a)
%
% See also: importdata, load, iLoad_ini, Loaders
%
% Part of: Loaders utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. 
% Version: $Date$
% (c) E.Farhi, ILL. License: EUPL.

% calls:    urlread
% optional: uigetfiles, looktxt, unzip, untar, gunzip (can do without)
% private:  iLoad_loader_auto, iLoad_config_load, iLoad_config_save, 
%           iLoad_import, iLoad_loader_check, iLoad_findfield

  persistent config

  if isempty(config)
    config  = iLoad_config_load;
    warning off backtrace % limit level of output messages 
  end

  data = []; format = [];
  if nargin == 0, filename=''; end
  if nargin < 2,  loader = ''; end
  if nargin ==1 && (ischar(filename) || isstruct(filename))
    if any(strcmp(filename, {'load config','config','force','force load config','formats','display config','load','save','compile','check'}))
      [data, format] = iLoad('', filename, varargin{:});
      return
    elseif  isstruct(filename)
      config = filename;
      loader = 'save';
      read_anytext('config');
    end
  end
  % is the loader in fact an option ?
  if ischar(loader) && (strncmp(loader,'--',2) || strncmp(loader,'-',1))
    varargin{end+1} = loader;
    loader='';
  end

  % ------------------------------------------------------------------------------
  % begin configuration stuff
  % ------------------------------------------------------------------------------

  if any(strcmp(loader, {'load config','config','load'}))
    
    % look for a specific importer when filename is specified
    if ~isempty(filename)
      data = {};
      for index=1:length(config.loaders)
        this = config.loaders{index};
        keepme=0;
        if ~isempty(strfind(lower(this.name), lower(filename))) | ~isempty(strfind(lower(this.method), lower(filename)))
          keepme = 1;
        elseif isfield(this,'extension') && any(strncmpi(filename, this.extension, length(filename)))
          keepme=1;
        end
        if keepme, data = { data{:} this }; end
      end
      if length(data) == 1
        data = data{1}; 
      end
    else
      data = config;
    end
    return
  elseif any(strcmp(loader, {'force','force load config','compile','check'}))
  
    config  = iLoad_config_load;                          % force omport
    data = iLoad_check_compile(filename, loader, config); % check executables
    return
  elseif strcmp(loader, 'formats') || strcmp(loader, 'display config')  || strcmp(loader, 'list')
    if ~isempty(filename)
      data    = iLoad(filename, 'load config', varargin{:});
    else
      data = iLoad('','load config', varargin{:});
    end
    data = iLoad_display_formats(data, config);
    return
  elseif any(strcmp(loader, {'save config','save'}))
    if isempty(filename) || nargin == 1
      config  = iLoad('','load config');
    else
      config = filename;
    end
    data = iLoad_config_save(config);
    return
  end
  % ------------------------------------------------------------------------------
  % end configuration stuff
  % ------------------------------------------------------------------------------

  % multiple file handling (directories, etc)
  if iscellstr(filename) & length(filename) > 1 & ~isempty(filename)
    data  = {};
    format= data;
    for index=1:numel(filename)
      [this_data, this_format] = iLoad(filename{index}, loader, varargin{:});
      if ~iscell(this_data),   this_data  ={ this_data }; end
      if ~iscell(this_format), this_format={ this_format }; end
      data  = [ data ; this_data(:) ];
      format= [ format ; this_format(:) ];
      % clear this_data this_format can not clear in for
    end
    return
  end

  if iscellstr(filename) & length(filename) == 1
    filename = filename{1};
  end

  url = false; % flag indicating that 'filename' is a temp file to be removed afterwards

  % handle single file name
  if ischar(filename) & length(filename) > 0
    % handle ~ substitution for $HOME
    if filename(1) == '~' && (length(filename==1) || filename(2) == '/' || filename(2) == '\')
      filename(1) = '';
      if usejava('jvm')
        filename = [ char(java.lang.System.getProperty('user.home')) filename ];
      elseif ~ispc  % does not work under Windows
        filename = [ getenv('HOME') filename ];
      end
    end

    % local/distant file (general case)
    f=find(filename == '#');
    if length(f) == 1 && f > 1  % the filename contains an internal link (HTML anchor)
      [fileroot,filesub]=strtok(filename, '#');
      if length(find(fileroot == '*')) || isdir(fileroot) % directory with token to search
        % get full recursive listing
        filerec=getAllFiles(fileroot);
        % find items that match the token
        found = ~cellfun('isempty',strfind(filerec, filesub(2:end)));
        filerec = filerec(found);
        if isempty(filerec), return; end
        [data, format]=iLoad(filerec, loader, varargin{:});
      else
        [data, format]=iLoad(fileroot, loader, varargin{:});
      end
      
      % now search pattern in the file names or fields
      if iscell(data) && (length(find(fileroot == '*')) || isdir(fileroot))
        this_data = {}; this_format = {};
        % name is a directory. We search for the pattern in file names
        f = [];
        for index=1:length(data)
          if any(strfind(data{index}.Source, filesub(2:end)))
            this_data{end+1}   = data{index};
            this_format{end+1} = loader;
          end
        end
        data   = this_data;
        format = this_format;
        return
      else
        f = iLoad_findfield(data, filesub(2:end)); % calls private function
        if iscell(f) && all(cellfun(@isempty, f)), f=[]; end
        if ~isempty(f)
          this_data = {}; this_format = {};
          % 'data' will be a cell array of structure based on the initial one
          for index=1:length(f)
            ret = data;
            fields = textscan(f{index},'%s','Delimiter','.'); % split the path in the structure with '.' char
            ret.Data = getfield(data,fields{1}{:});             % access that path by expanding it;
            try
              fields{1}{1} = 'Headers';
              ret.Headers = getfield(data,fields{1}{:});
            end
            this_data{end+1}   = ret;
            this_format{end+1} = loader;

          end   
          data   = this_data;
          format = this_format;
          return
        else
          fprintf(1, 'iLoad: Warning: Could not find pattern "%s". Importing whole file...\n', filesub(2:end));
          filename = fileroot;
        end
      end
    end
    
    if strncmp(filename, 'file://', length('file://'))
      filename = filename(7:end); % remove 'file://' from local name
    end
    
    % handle / to \ substitution for Windows systems, not in URLs
    if ~(strncmp(filename, 'http://', length('http://')) | ...
         strncmp(filename, 'https://',length('https://'))   | ...
         strncmp(filename, 'ftp://',  length('ftp://'))   | ...
         strncmp(filename, 'file://', length('file://')) )
      if    ~ispc, filename = strrep(filename, '\', filesep);
      elseif ispc, filename = strrep(filename, '/', filesep);
      end
    end
    
    if isempty(dir(filename)) && ~isempty(dir(fullfile(ifitpath,'Data',filename)))
          filename = fullfile(ifitpath,'Data',filename);
    end
    
    if isdir(filename), filename = [ filename filesep '*']; end % all elements in case of directory
    
    % handle single file name (possibibly with wildcard)
    if ~isempty(find(filename == '*')) | ~isempty(find(filename == '?'))  % wildchar !!#
      [filepath,name,ext]=fileparts(filename);  % 'file' to search
      if isempty(filepath), filepath = pwd; end
      this_dir = dir(filename);
      if isempty(this_dir), return; end % directory is empty
      % removes '.' and '..'
      index    = find(~strcmp('.', {this_dir.name}) & ~strcmp('..', {this_dir.name}));
      this_dir = char(this_dir.name);
      this_dir = (this_dir(index,:));
      if isempty(this_dir), return; end % directory only contains '.' and '..'
      rdir     = cellstr(this_dir); % original directory listing as cell
      rdir     = strcat([ filepath filesep ], char(rdir));
      filename = cellstr(rdir);
      [data, format] = iLoad(filename, loader, varargin{:});
      return
    end
    
    % handle file on the internet
    if strncmp(filename, 'http://', length('http://')) ...
     | strncmp(filename, 'https://',length('https://')) ...
     | strncmp(filename, 'ftp://',  length('ftp://'))
      tmpfile = tempname;
      % Keep file extension, may be useful for iData load
      [filepath,name,ext] = fileparts(filename);
      tmpfile = [tmpfile ext];
      use_wget = false;
      if ~usejava('jvm')
        use_wget = true;
      else
        % access the net. Proxy settings must be set (if any).
        try
          % write to temporary file
          tmpfile = urlwrite(filename, tmpfile);
        catch ME
          disp(ME.message);
          use_wget = true;
        end
      end
      if use_wget
        % Fall back to using wget
        cmd = ['wget ' filename ' -O ' tmpfile]; disp(cmd)
        [status, result] = system(cmd);
        if status
          disp(result);
          error([ mfilename ': Can not get URL ' filename ]);
        end
      end
      filename = tmpfile;
    end
    
    % handle compressed files (local or distant)
    [pathstr, name, ext] = fileparts(filename);
    if     strcmp(ext, '.zip'), cmd = 'unzip';
    elseif strcmp(ext, '.tar'), cmd = 'untar';
    elseif strcmp(ext, '.gz') || strcmp(ext, '.tgz'),  cmd = 'gunzip';
    elseif strcmp(ext, '.Z'),   cmd = 'uncompress';
    else                        cmd=''; end
    if ~isempty(cmd)
      % this is a compressed file/url. Extract to temporary dir.
      if strcmp(cmd, 'uncompress')
        copyfile(filename, tempdir, 'f');
        try
          system(['uncompress ' tempdir filesep name ext ]);
          filename = [ tempdir filesep name ];
          url = true;
        catch ME
          fprintf(1, 'iLoad: Can''t extract file "%s" (%s).\n', filename,cmd);
          warning(ME.message)
          return
        end
      elseif exist(cmd)
        % extract to temporary dir
        try
          filenames = feval(cmd, filename, tempdir);
        catch ME
          fprintf(1, 'iLoad: Can''t extract file "%s" (%s).\n', filename,cmd);
          warning(ME.message)
          return
        end
        [data, format] = iLoad(filenames, loader, varargin{:}); % is now local
        for index=1:length(filenames)
          try
            delete(filenames{index});
          catch ME
            fprintf(1,'iLoad: Can''t delete temporary file "%s" (%s).\n', filename{index},cmd);
            warning(ME.message)
          end
        end
        return
      end
    end
    
    % The import takes place HERE ================================================
    if isdir(filename), filename = [ filename filesep '*']; end % all elements in case of directory
    % handle the '%20' character replacement as space
    filename = strrep(filename, '%20',' ');
    while filename(end) == ';', filename(end)=''; end % in case there is a leading ';' in place of \n
    try
      [data, format] = iLoad_import(filename, loader, varargin{:}); % <<<<<<<<<<
    catch ME
      fprintf(1, 'iLoad: Failed to import file %s. Ignoring.\n  %s\n', filename, ME.message);
      data=[];
    end
  elseif isempty(filename)
    config = iLoad('','load config');
    % build the list of supported data formats
    filterspec = {'*.*', 'All files (*.*)' };
    config = iLoad('','load config'); % persistent
    for index=1:numel(config.loaders)
      this = config.loaders{index};
      if isfield(this, 'name') && isfield(this,'method')
        if isfield(this, 'extension') && ~isempty(this.extension)
          ext=this.extension;
        else ext='*'; end
        if ischar(ext)  ext = [ '*.' ext ]; end
        if iscell(ext), ext = sprintf('*.%s;',ext{:}); end
        filterspec(end+1,1:2) = { ext, [ this.name ' (' ext ')' ]};
      end
    end
    if isempty(filterspec), filterspec='*.*'; end
    if exist('uigetfiles') && ((isfield(config, 'UseSystemDialogs') && strcmp(config.UseSystemDialogs, 'no')) || ~usejava('jvm'))
        [filename, pathname] = uigetfiles(filterspec,'Select file(s) to load');
    else
      if usejava('swing')
        setappdata(0,'UseNativeSystemDialogs',false);
        [filename, pathname] = uigetfile(filterspec, 'Select file(s) to load', 'MultiSelect', 'on');
      else
        [filename, pathname] = uigetfile(filterspec, 'Select a file to load');
      end
    end
    if isempty(filename),    return; end
    if isequal(filename, 0), return; end
    filename = strcat(pathname, filesep, filename);
    if ~iscellstr(filename)
      if isdir(filename)
        filename = [ filename filesep '*']; 
      end % all elements in case of directory
    end
    [data, format] = iLoad(filename, loader, varargin{:});
  else
    % data not empty, but not a file name
    data = iLoad_loader_check([ inputname(1) ' variable of class ' class(filename) ], filename, 'variable');
    format= '' ;
  end

  % remove temporary file if needed
  if (url)
    try
    delete(filename);
    catch ME
    fprintf(1,'iLoad: Can''t delete temporary file "%s".\n', filename);
    warning(ME.message);
    end
  end
  
  % NESTED functions so that we can access 'config' faster
  % ----------------------------------------------------------------------------
  % private/nested function to import single data with given method(s)
  function [data, loader] = iLoad_import(filename, loader, varargin)
    data = []; isbinary=0;
    verbose = 1; % set this to 1 to get more output for debugging
    
    if isempty(dir(filename))
      loader = 'Failed to load file (does not exist)';
      return
    end
    if isempty(loader), loader='auto'; end
    if strcmp(loader, 'auto')
      [loader, isbinary] = iLoad_loader_auto(filename);
    elseif strcmp(loader, 'gui')
      [dummy, filename_short, ext] = fileparts(filename);
      filename_short = [ filename_short ext];
      [loader, isbinary] = iLoad_loader_auto(filename);
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
          data = iLoad_import(filename, this_loader, varargin{:});
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
    varg = {};
    if iscell(loader.options)
      varg = { filename, loader.options{:}, varargin{:} };
    elseif ischar(loader.options)
      if ~isempty(loader.options), loader.options = [ ' ' loader.options ]; end
      varg = { [ filename loader.options ], varargin{:} };
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
  
  % -----------------------------------------------------------
  % private/nested function to determine which parser to use to analyze content
  function [loaders, isbinary] = iLoad_loader_auto(file)
  
    verbose = 0;  % set it to 1 to display Loader auto-selection from extension/pattern

    % get the configuration
    % config  = iLoad('','load config');
    
    loaders = config.loaders;
    isbinary= 0;
    % read start of file
    [fid,message] = fopen(file, 'r');
    if fid == -1
      fprintf(1, 'iLoad: %s: %s. Check existence/permissions.\n', file, message );
      error([ 'Could not open file ' file ' for reading. ' message '. Check existence/permissions.' ]);
    end
    file_start = fread(fid, 1000, 'uint8=>char')';
    
    % loop to test each format for patterns
    formats = loaders;
    loaders ={}; % will hold those formats which are applicable
    loaders_count=0;
    % identify by extensions
    [dummy, dummy, fext] = fileparts(file);
    if ~isempty(fext) && fext(1)=='.', fext(1)=[]; end
    % check if this is a text file
    if length(find(file_start >= 32 & file_start < 127))/length(file_start) < 0.4
      isbinary = 1; % less than 90% of printable characters
    else
      % clean file start with spaces, remove EOL
      file_start = [ file_start fread(fid, 9000, 'uint8=>char')' ];
      file_start(isspace(file_start)) = ' ';
    end
    fclose(fid);
    
    % identify by patterns
    for index=1:length(formats)
      loader = formats{index};
      if ~isstruct(loader), continue; end

      if ~isfield(loader,'patterns'), loader.patterns=nan; end
      
      patterns_found = 0;
      if verbose, fprintf(1,'iLoad: method %s: file %s: analyzing\n', loader.name, file); end
      
      % the loader is selected if: 
      %   loader has extension and extension matches, no patterns to search
      if ~patterns_found && (isempty(loader.patterns) || isnumeric(loader.patterns))
        if ~isfield(loader,'extension'), ext=''; 
        else ext=loader.extension; end
        if ischar(ext) && ~isempty(ext), ext={ ext }; end
        if ~isempty(ext) && ~isempty(fext) 
          if any(strcmpi(fext, ext))
            patterns_found  = 1;  % extension does match
            if verbose, disp([ 'iLoad: method ' loader.name ': ' file ': extension ' ext{1} ' matches' ]); end
          end
        end
      end
        
      %   loader has patterns and not binary and patterns match file_start, 
      %   whatever be the extension
      % ~isbinary may be removed in case it suppresses e.g. text/binary formats such as EDF, ADSC, ...
      if ~isempty(loader.patterns) && ischar(loader.patterns)
        loader.patterns=cellstr(loader.patterns); 
      end
      
      if ~patterns_found && ~isempty(loader.patterns) && iscellstr(loader.patterns) && ~isbinary 
      
        % check all patterns in text file
        all_match = 1;
        for index_pat=1:numel(loader.patterns) % all patterns must match
          if isempty(regexpi(file_start, loader.patterns{index_pat}, 'once'))
            all_match=0;     % at least one pattern does not match
            if verbose, fprintf(1,'iLoad: method %s: file %s: at least one pattern does not match (%s)\n', loader.name, file, loader.patterns{index_pat}); end
            break;
          end
        end % for patterns
        if all_match, patterns_found=1; 
          if verbose, disp([ 'iLoad: method ' loader.name ': ' file ': patterns match' ]); end
        end
      end
        
      %   loader has no extension and no patterns
      if ~patterns_found && isempty(ext) && isempty(loader.patterns)
        patterns_found = 1; % we will try all non selective loaders
        if verbose, disp([ 'iLoad: method ' loader.name ': ' file ': no patterns nor extension: select by default' ]); end
      end
        
      if patterns_found
        loaders_count = loaders_count+1;
        loaders{loaders_count} = loader;
      end

    end % for index

  end % iLoad_loader_auto

end % iLoad (main)


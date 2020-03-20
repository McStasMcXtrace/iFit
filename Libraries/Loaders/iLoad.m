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
%   iLoad([1 2 3 4])
%   iLoad(structure)
%   iLoad({ 1:10, 3:4:80 })
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
%   iLoad verbose       switch to verbose mode (more messages)
%   iLoad silent        switch to silent mode (no messages)
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
% Author:  E. Farhi <emmanuel.farhi@synchrotron-soleil.fr>. 
% Version: $Date$ $Version$ $Author$
% 

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
    if any(strcmp(filename, {'load config','config','force','force load config','formats','display config','load','save','compile','check','silent','verbose'}))
      [data, format] = iLoad('', filename, varargin{:});
      return
    elseif  isconfig(filename)
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
  
    config  = iLoad_config_load;                          % force import
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
  elseif any(strcmp(loader, {'silent','verbose'}))
    if strcmp(loader, 'silent')
      config.verbosity = 0;
    elseif strcmp(loader, 'verbose')
      config.verbosity = 2;
    end
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

    % local/distant file (general case) #token 
    f=find(filename == '#');
    if length(f) == 1 && f > 1  % the filename contains an internal link (HTML anchor)
      [fileroot,filesub]=strtok(filename, '#');
      if length(find(fileroot == '*')) || isdir(fileroot) % directory with #token to search
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
      
      % now search wildcard in the file names or fields
      if iscell(data) && (~isempty(find(fileroot == '*')) || isdir(fileroot))
        this_data = {}; this_format = {};
        % name is a directory. We search for the pattern in file names
        f = [];
        for index=1:length(data)
          if any(strfind(data{index}.Source, filesub(2:end))) % #token
            this_data{end+1}   = data{index};
            this_format{end+1} = loader;
          end
        end
        data   = this_data;
        format = this_format;
        return
      else
        % search #token in imported data
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
    end % filename with '#' anchor
    
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
    if isempty(dir(filename))
      warning([ mfilename ': possibly invalid filename ' filename ]);
    end
    % handle the '%20' character replacement as space
    filename = strrep(filename, '%20',' ');
    while filename(end) == ';', filename(end)=''; end % in case there is a leading ';' in place of \n
    try
      [data, format] = iLoad_import(filename, loader, config, varargin{:}); % <<<<<<<<<<
    catch ME
      disp(getReport(ME))
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
  elseif isiloadstruct(filename)
    data   = filename;
    format = filename.Loader;
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

end % iLoad (main)

% ------------------------------------------------------------------------------
function tf = isconfig(filename)
  % ISCONFIG Return true when filename is a struct with an iLoad config
  tf = false;
  if ~isstruct(filename), return; end
  for f = {'FileName','UseSystemDialogs','MeX','loaders'}
    if ~isfield(filename, f{1}), return; end
  end
  if ~iscell(filename.loader), return; end
  tf = true;
end

function tf = isiloadstruct(filename)
  % ISILOADSTRUCT Return true when this is already an iLoad struct
  
  tf = false;
  for f={'Data' 'Source' 'Title' 'Date' 'Format' 'Command' 'Creator' 'User' 'Label' 'Loader'}
    if ~isfield(filename, f{1}), return; end
  end
  tf = true;
end


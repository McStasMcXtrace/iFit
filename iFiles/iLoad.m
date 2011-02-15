function [data, format] = iLoad(filename, loader)
% [data, loader] = iLoad(file, loader)
%
% imports any data into Matlab. 
% The definition of specific formats can be set in the iLoad_ini.m file.
% These formats can be obtained using [config, configfile]=iLoad('','load config').
% the iLoad_ini configuration file can be saved in the Preference directory
% using [config, configfile] = iLoad(config,'save config').
% A list of all supported formats is shown with iLoad('formats');
%
%   Default supported formats include: any text based including CSV, Lotus1-2-3, SUN sound, 
%     WAV sound, AVI movie, NetCDF, FITS, XLS, BMP GIF JPEG TIFF PNG ICO images,
%     HDF4, HDF5, MAT workspace, XML.
%   Other specialized formats include: McStas, ILL, SPEC, ISIS/SPE, INX, EDF.
%   Compressed files are also supported, with on-the-fly extraction (zip, gz, tar, Z).
%   Distant files are supported through e.g. URLs such as 
%     file://, ftp:// and http://
%
% input arguments:
%   file:   file name, or cell of file names, or any Matlab variable, or a URL
%             or an empty string (then pops'up a file selector)
%   loader: a function name to use as import routine, OR a structure with:
%             loader = 'auto' (default) 
%             loader = 'gui' (asks for the format to use)
%             loader.method = 'function name'
%             loader.options= string of options to catenate after file name
%                          OR cell of options to use as additional arguments
%                             to the method
%
% output variables:
%   data:   a single structure containing file data, or a cell of structures
%   loader: the loader that was used for importation, or a cell of loaders.
%
% example: iLoad; iLoad('file');
%
% See also: importdata, load, iLoad_ini
%
% Part of: iFiles utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. % Version: $Revision: 1.34 $

% calls:    urlread
% optional: uigetfiles, looktxt, unzip, untar, gunzip (can do without)
% private:  iLoad_loader_auto, iLoad_config_load, iLoad_config_save, iLoad_import, iLoad_loader_check

persistent config

data = []; format = [];
if nargin == 0, filename=''; end
if nargin < 2,  loader = ''; end
if strcmp(loader, 'load config') || strcmp(filename, 'load config')
  if isempty(config), config  = iLoad_config_load; end
  data = config;
  return
elseif strcmp(loader, 'force load config') || strcmp(filename, 'force load config')
  config  = iLoad_config_load;
  data = config;
  return
elseif strcmp(loader, 'formats') || strcmp(filename,'formats') || strcmp(loader, 'display config')
  data = iLoad('','load config');
  fprintf(1, ' EXT                    READER  DESCRIPTION\n');
  fprintf(1, '-----------------------------------------------------------------\n');  
  for index=1:length(data.loaders)
    this=data.loaders{index};
    if isfield(this,'postprocess'), 
      if ~isempty(this.postprocess), this.method = [ this.method '/' this.postprocess ]; end
    end
    if length(this.method)>25, this.method = [ this.method(1:22) '...' ]; end
    if ~isfield(this,'extension'), this.extension = '*';
    elseif isempty(this.extension), this.extension='*'; end
    if iscellstr(this.extension)
      fprintf(1,'%4s %25s  %s\n', upper(this.extension{1}), this.method,this.name);
      for j=2:length(this.extension),fprintf(1,'  |.%s\n', upper(this.extension{j})); end
    else
      fprintf(1,'%4s %25s  %s\n', upper(this.extension), this.method,this.name);
    end
  end
  disp([ '% iLoad configuration file: ' config.FileName ]);
  return
elseif strcmp(loader, 'save config') || strcmp(filename, 'save config')
  if isempty(filename)
    config  = iLoad('','load config');
  else
    config = filename;
  end
  data = iLoad_config_save(config);
  return
end

% multiple file handling
if iscellstr(filename) & length(filename) > 1 & ~isempty(filename)
  data  = cell(length(filename(:)), 1);
  format= data;
  for index=1:length(filename(:))
    [data{index}, format{index}] = iLoad(filename{index}, loader);
  end
  return
end

if iscellstr(filename) & length(filename) == 1
  filename = filename{1};
end

% handle single file name
if ischar(filename) & length(filename) > 0
  % local file (general case)
  if strncmp(filename, 'file://', length('file://'))
    filename = filename(8:end); % remove 'file://' from name
  end
  
  % handle / to \ substitution for Windows systems, not in URLs
  if ~(strncmp(filename, 'http://', length('http://')) | ...
       strncmp(filename, 'ftp://', length('ftp://'))   | ...
       strncmp(filename, 'file://', length('file://')) )
    if    ~ispc, filename = strrep(filename, '\', filesep);
    elseif ispc, filename = strrep(filename, '/', filesep);
    end
  end
  
  if isdir(filename), filename = [ filename filesep '*']; end % all elements in case of directory
  
  % handle single file name (possibibly with wildcard)
  if ~isempty(find(filename == '*')) | ~isempty(find(filename == '?'))  % wildchar !!#
    [filepath,name,ext,versn]=fileparts(filename);  % 'file' to search
    if isempty(filepath), filepath = pwd; end
    this_dir = dir(filename);
    if isempty(this_dir), return; this_dir = dir(filepath); end
    index = find(real([this_dir.isdir]) == 0);
    this_dir = char(this_dir.name);
    this_dir = (this_dir(index,:));
    rdir = cellstr(this_dir); % original directory listing as cell
    rdir = strcat([ filepath filesep ], char(rdir));
    filename = cellstr(rdir);
    [data, format] = iLoad(filename, loader);
    return
  end
  
  % handle compressed files
  [pathstr, name, ext] = fileparts(filename);
  if     strcmp(ext, '.zip'), cmd = 'unzip';
  elseif strcmp(ext, '.tar'), cmd = 'untar';
  elseif strcmp(ext, '.gz'),  cmd = 'gunzip';
  elseif strcmp(ext, '.Z'),   cmd = 'uncompress';
  else                       cmd=''; end
  if ~isempty(cmd)
    % this is a compressed file/url. Extract to temporary dir.
    if strcmp(cmd, 'uncompress')
      copyfile(filename, tempdir, 'f');
      try
        system(['uncompress ' tempdir filesep name ext ]);
        [data, format] = iLoad([ tempdir filesep name ], loader); % is now local in tempdir
      catch
      end
      delete([ tempdir filesep name ]);
      return
    elseif exist(cmd)
      % extract to temporary dir
      filenames = feval(cmd, filename, tempdir);
      [data, format] = iLoad(filenames, loader); % is now local
      for index=1:length(filenames)
        delete(filenames{index});
      end
      return
    end
  end
  
  % handle files on the internet
  if strncmp(filename, 'http://', length('http://')) | strncmp(filename, 'ftp://', length('ftp://'))
    % access the net. Proxy settings must be set (if any).
    data = urlread(filename);
    % write to temporary file
    tmp_file = tempname;
    fid = fopen(tmp_file, 'w');
    fprintf(fid, '%s', data);
    fclose(fid);
    [data, format] = iLoad(tmp_file, loader); % is now local
    delete(tmp_file);
    return
  else
    % local file (general case)
    if strncmp(filename, 'file://', length('file://'))
      filename = filename(8:end); % remove 'file://' from name
    end
    % The import takes place HERE
    if isdir(filename), filename = [ filename filesep '*']; end % all elements in case of directory
    [data, format] = iLoad_import(filename, loader);
  end
elseif isempty(filename)
  config = iLoad('','load config');
  if exist('uigetfiles') & strcmp(config.UseSystemDialogs, 'no')
      [filename, pathname] = uigetfiles('.*','Select file(s) to load');
  else
    if usejava('swing')
      setappdata(0,'UseNativeSystemDialogs',false);
      [filename, pathname] = uigetfile('*.*', 'Select file(s) to load', 'MultiSelect', 'on');
    else
      [filename, pathname] = uigetfile('*.*', 'Select a file to load');
    end
  end
  if isempty(filename), return; end
  filename = strcat(pathname, filesep, filename);
  if ~iscellstr(filename)
    if isdir(filename)
      filename = [ filename filesep '*']; 
    end % all elements in case of directory
  end
  [data, format] = iLoad(filename, loader);
else
  % data not empty, but not a file name
  data = iLoad_loader_check([ inputname(1) ' variable of class ' class(filename) ], filename, 'variable');
  format= '' ;
end

% -----------------------------------------------------------
% private function to import single data with given method(s)
function [data, loader] = iLoad_import(filename, loader)
  data = [];
  if isempty(loader), loader='auto'; end
  if strcmp(loader, 'auto')
    loader = iLoad_loader_auto(filename);
  elseif strcmp(loader, 'gui')
    [dummy, filename_short, ext] = fileparts(filename);
    filename_short = [ filename_short ext];
    loader = iLoad_loader_auto(filename);
    loader_names=[loader{:}];
    tmp         = cell(size(loader_names)); tmp(:)={' - '};
    loader_names= strcat({loader_names.name}, tmp,{loader_names.method});
    loader_index=listdlg(...
      'PromptString',...
        {'Select suitable import methods',['to import file ' filename_short ]}, ...
      'SelectionMode','Multiple',...
      'ListString', loader_names, ...
      'ListSize', [300 160], ...
      'Name', ['Loader for ' filename_short ]);
    if isempty(loader_index), loader=[]; return; end
    loader=loader(loader_index);
  elseif ischar(loader)
    % test if loader is the user name of a function
    config = iLoad('','load config');
    formats = config.loaders;
    loaders={};
    loaders_count=0;
    for index=1:length(formats)
      this_loader = formats{index};
      if ~isempty(strfind(this_loader.name, loader)) || ~isempty(strfind(this_loader.method, loader))
        loaders_count = loaders_count+1;
        loaders{loaders_count} = this_loader;
      end
    end
    if ~isempty(loaders) loader = loaders; end
  end
  
  % handle multiple loaders (cell or struct array)
  if (iscell(loader) | isstruct(loader)) & length(loader) > 1
    loader=loader(:);
    for index=1:length(loader)
      if iscell(loader), this_loader = loader{index};
      else this_loader = loader(index); end
      try
        data = iLoad_import(filename, this_loader);
      catch
        disp(lasterr)
        [dummy, name_short, ext] = fileparts(filename);
        fprintf(1, 'iLoad: Failed to import file %s with method %s (%s). Ignoring.\n', name_short, this_loader.name, this_loader.method);
        data = [];
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

  % handle single char loaders (IMPORT takes place HERE)
  if ischar(loader)
    tmp=loader; clear loader;
    loader.method = tmp; loader.name=tmp; loader.options='';
  end
  if isempty(loader.method), return; end
  fprintf(1, 'iLoad: Importing file %s with method %s (%s)\n', filename, loader.name, loader.method);
  if isempty(loader.options)
    data = feval(loader.method, filename);
  elseif iscell(loader.options)
    data = feval(loader.method, filename, loader.options{:})
  elseif ischar(loader.options)
    try
    data = feval(loader.method, [ filename ' '  loader.options ]);
    catch
    data = feval(loader.method, filename,loader.options);
    end
  end
  data = iLoad_loader_check(filename, data, loader);
  if isfield(loader, 'name') data.Format = loader.name; 
  else data.Format=[ loader.method ' import' ]; end
  return
  
% -----------------------------------------------------------
% private function to make the data pretty looking
function data = iLoad_loader_check(file, data, loader)

  % handle case when a single file generates a data set
  if isstruct(data) & length(data)>1
    for index=1:length(data)
      data(index) = iLoad_loader_check(file, data(index), loader);
    end
    return
  elseif iscellstr(data)
    fprintf(1, 'iLoad: Failed to import file %s with method %s (%s). Got a cell of strings. Ignoring\n', file, loader.name, loader.method);
  elseif iscell(data) & length(data)>1
    newdata=[];
    for index=1:length(data)
      newdata(index) = iLoad_loader_check(file, data{index}, loader);
    end
    data = newdata; % now an array of struct
    return
  end
  
  name='';
  if isstruct(loader),
    method = loader.method;
    options= loader.options;
    if isfield(loader, 'name'), name=loader.name; end
  else
    method = loader; options=''; 
  end

  if isempty(method), method='iData/load'; end
  if strcmp(loader, 'variable')
    method='iData/load';
  end
  if isempty(name), name=method; end
  if iscell(options), options= cellstr(options{1}); options= [ options{1} ' ...' ]; end
  if ~isfield(data, 'Source')  & ~isfield(data, 'Date') & ~isfield(data, 'Format') ...
   & ~isfield(data, 'Command') & ~isfield(data,' Data')
    new_data.Data = data;
    % transfer some standard fields as possible
    if isfield(data, 'Source'), new_data.Source = data.Source; end
    if isfield(data, 'Title'), new_data.Title = data.Title; end
    if isfield(data, 'Date'), new_data.Date = data.Date; end
    if isfield(data, 'Label'), new_data.Label = data.Label; end
    
    data = new_data;
    
  end
  
  if ~isfield(data, 'Source') && ~isfield(data, 'Filename'),  data.Source = file;
  elseif isfield(data, 'Filename'), data.Source = data.Filename; end

  if ~isfield(data, 'Title'),   
    [pathname, filename, ext] = fileparts(file);
    if ~strcmp(loader, 'variable'), data.Title  = [ 'File ' filename ext ' ' name  ];
    else data.Title  = [ 'File ' filename ext ]; end
  end
  
  if ~isfield(data, 'Date')
    if strcmp(loader, 'variable') data.Date   = datestr(now); 
    else d=dir(file); data.Date=d.date; end
  end

  if ~isfield(data, 'Format'),
    if ~strcmp(loader, 'variable'), data.Format  = [ name ' import with Matlab ' method ];  
    else data.Format  = [ 'Matlab ' method ]; end
  end
  if ~isfield(data, 'Command'),
    if strcmp(loader, 'variable')
      data.Command = [ method '(' file ', '''  options ''')' ];
    else
      data.Command = [ method '(''' file ''', '''  options ''')' ];
    end
  end
  if ~isfield(data, 'Creator'), data.Creator = [ name ' iData/load/' method ]; 
  else data.Creator = [ name ' iData/load/' method ' - ' data.Creator ]; end
  if ~isfield(data, 'User'),
    if isunix
      data.User    = [ getenv('USER') ' running on ' computer ' from ' pwd ];
    else
      data.User    = [ 'User running on ' computer ' from ' pwd ];
    end
  end
  return

% -----------------------------------------------------------
% private function to determine which parser to use to analyze content
% if allformats == 1, no pattern search is done
function loaders = iLoad_loader_auto(file)
  config  = iLoad('','load config');
  loaders = config.loaders;
    
  % read start of file
  [fid,message] = fopen(file, 'r');
  if fid == -1
    error([ 'Could not open file ' file ' for reading. ' message '. Check existence/permissions.' ]);
  end
  file_start = fread(fid, 10000, 'uint8=>char')';
  fclose(fid);
  % loop to test each format for patterns
  formats = loaders;
  loaders={};
  loaders_count=0;
  % identify by extensions
  [dummy, dummy, fext] = fileparts(file);
  fext=strrep(fext,'.','');
  
  % identify by patterns
  for index=1:length(formats)
    loader = formats{index};
    if ~isstruct(loader), break; end

    if exist(loader.method)
      if strcmp(loader.method, 'looktxt') && ...
              length(find(file_start >= 32 & file_start < 127))/length(file_start) < 0.9
        % fprintf(1,'iLoad: skip method %s as file %s is probably binary\n', loader.method, file);
        patterns_found  = 0;
        continue;  % does not use looktxt for binary data files
      end
      patterns_found  = 1;
      if ~isfield(loader,'patterns') loader.patterns=''; end
      if isempty(loader.patterns)  % no pattern to search, test extension
        if ~isfield(loader,'extension'), ext=''; 
        else ext=loader.extension; end
        if ischar(ext) && length(ext), ext=cellstr(ext); end
        if length(ext) && length(fext) 
          if isempty(strmatch(lower(fext), lower(ext), 'exact'))
            patterns_found  = 0;  % extension does not match
            % fprintf(1,'iLoad: method %s file %s: extension does not match (%s) \n', loader.name, file, fext);
          end
        else
          patterns_found  = 1;    % no extension, no patterns: try loader anyway
        end
      else  % check patterns
        if ischar(loader.patterns), loader.patterns=cellstr(loader.patterns); end
        for index_pat=1:length(loader.patterns(:))
          if isempty(strfind(file_start, loader.patterns{index_pat}))
            patterns_found=0;     % at least one pattern does not match
            % fprintf(1,'iLoad: method %s file %s: at least one pattern does not match (%s)\n', loader.name, file, loader.patterns{index_pat});
            continue;
          end
        end % for patterns
      end % if patterns
      
      if patterns_found
        loaders_count = loaders_count+1;
        loaders{loaders_count} = loader;
      end
    end % if exist(method)
  end % for index
  
  return;

% -----------------------------------------------------------
% private function to save the configuration and format customization
function config = iLoad_config_save(config)
  data = config.loaders;
  format_names  ={};
  format_methods={};
  format_unique =ones(1,length(data));
  % remove duplicated format definitions
  for index=1:length(data)
    if ~isempty(data{index}.name) & ~isempty(strmatch(data{index}.name, format_names, 'exact')) & ~isempty(strmatch(data{index}.method, format_methods, 'exact'))
      format_unique(index) = 0; % already exists. Skip it.
      format_names{index}  = '';
      format_methods{index}= '';
    else
      format_names{index} = data{index}.name;
      format_methods{index} = data{index}.method;
    end
  end
  data = data(find(format_unique));
  config.loaders = data;
  % save iLoad.ini configuration file
  % make header for iLoad.ini
  config.FileName=fullfile(prefdir, 'iLoad.ini'); % store preferences in PrefDir (Matlab)
  str = [ '% iLoad configuration script file ' sprintf('\n') ...
          '%' sprintf('\n') ...
          '% Matlab ' version ' m-file ' config.FileName sprintf('\n') ...
          '% generated automatically on ' datestr(now) ' with iLoad('''',''save config'');' sprintf('\n') ...
          '%' sprintf('\n') ...
          '% The configuration may be specified as:' sprintf('\n') ...
          '%     config = { format1 .. formatN }; (a single cell of format definitions, see below).' sprintf('\n') ...
          '%   OR a structure' sprintf('\n') ...
          '%     config.loaders = { format1 .. formatN }; (see below)' sprintf('\n') ...
          '%     config.UseSystemDialogs=''yes'' to use built-in Matlab file selector (uigetfile)' sprintf('\n') ...
          '%                             ''no''  to use iLoad file selector           (uigetfiles)' sprintf('\n') ...
          '%' sprintf('\n') ...
          '% User definitions of specific import formats to be used by iLoad' sprintf('\n') ...
          '% Each format is specified as a structure with the following fields' sprintf('\n') ...
          '%   method:   function name to use, called as method(filename, options...)' sprintf('\n') ...
          '%   extension:a single or a cellstr of extensions associated with the method' sprintf('\n') ...
          '%   patterns: list of strings to search in data file. If all found, then method' sprintf('\n') ...
          '%             is qualified' sprintf('\n') ...
          '%   name:     name of the method/format' sprintf('\n') ...
          '%   options:  additional options to pass to the method.' sprintf('\n') ...
          '%             If given as a string they are catenated with file name' sprintf('\n') ...
          '%             If given as a cell, they are given to the method as additional arguments' sprintf('\n') ...
          '%   postprocess: function called from iData/load after file import, to assign aliases, ...' sprintf('\n') ...
          '%             called as iData=postprocess(iData)' sprintf('\n') ...
          '%' sprintf('\n') ...
          '% all formats must be arranged in a cell, sorted from the most specific to the most general.' sprintf('\n') ...
          '% Formats will be tried one after the other, in the given order.' sprintf('\n') ...
          '% System wide loaders are tested after user definitions.' sprintf('\n') ...
          '%' sprintf('\n') ...
          '% NOTE: The resulting configuration must be named "config"' sprintf('\n') ...
          '%' sprintf('\n') ...
          class2str('config', config) ];
  [fid, message]=fopen(config.FileName,'w+');
  if fid == -1
    warning(['Error opening file ' config.FileName ' to save iLoad configuration.' ]);
    config.FileName = [];
  else
    fprintf(fid, '%s', str);
    fclose(fid);
    disp([ '% Saved iLoad configuration into ' config.FileName ]);
  end
  
% -----------------------------------------------------------
% private function to load the configuration and format customization
function config = iLoad_config_load
  loaders      = {};
  % read user list of loaders which is a cell of format descriptions
  if exist(fullfile(prefdir, 'iLoad.ini'), 'file')
    % there is an iLoad_ini in the Matlab preferences directory: read it
    configfile = fullfile(prefdir, 'iLoad.ini');
    fid = fopen(configfile, 'r');
    content = fread(fid, Inf, 'uint8=>char');
    fclose(fid);
    % evaluate content of file
    config=[]; eval(content(:)'); % this makes a 'config' variable
    if iscell(config)
      loaders = config; config=[];
      config.loaders = loaders;
      config.FileName= configfile;
    end
    disp([ '% Loaded iLoad format descriptions from ' config.FileName ]);
  elseif exist('iLoad_ini')
    config = iLoad_ini;
  end
  
  % check if other configuration fields are present, else defaults
  if ~isfield(config, 'UseSystemDialogs'), config.UseSystemDialogs = 'no'; end
  if ~isfield(config, 'FileName'),         config.FileName = ''; end
  
  loaders = config.loaders;
  
  % ADD default loaders: method, ext, name, options
  % default importers, when no user specification is given. 
  % These do not have any pattern recognition or postprocess
  formats = {...
    { 'looktxt', '',    'Data (text format with fastest import method)',    '--headers --binary --fast --comment= '}, ...
    { 'looktxt', '',    'Data (text format with fast import method)',       '--headers --binary --comment= '}, ...
    { 'looktxt', '',    'Data (text format)',                               '--headers --comment= '}, ...
    { 'wk1read', 'wk1', 'Lotus1-2-3 (first spreadsheet)',''}, ...
    { 'auread',  'au',  'NeXT/SUN (.au) sound',''}, ...
    { 'wavread', 'wav'  'Microsoft WAVE (.wav) sound',''}, ...
    { 'aviread', 'avi', 'Audio/Video Interleaved (AVI) ',''}, ...
    { 'mcdfread', {'nc','cdf'}, 'NetCDF 2 (.nc)',''}, ...
    { 'netcdf',   {'nc','cdf'}, 'NetCDF 1.0 (.nc)','','','load_netcdf1'}, ...
    { 'mfitsread',{'fits','fts'},'FITS',''}, ...
    { 'xlsread', 'xls', 'Microsoft Excel (first spreadsheet, .xls)',''}, ...
    { 'mimread',  {'bmp','jpg','jpeg','tiff','png','ico'}, 'Image/Picture',''}, ...
    { 'mhdf4read',{'hdf4','h4','hdf'},  'HDF4 raster image',''}, ...
    { 'hdf5extract',{'hdf5','h5'}, 'HDF5',''}, ...
    { 'load',    'mat', 'Matlab workspace (.mat)',''}, ...
    { 'csvread', 'csv', 'Comma Separated Values (.csv)',''}, ...
    { 'dlmread', 'dlm', 'Numerical single block',''}, ...
    { 'xmlread', 'xml', 'XML',''}, ...
    { 'importdata','',  'Matlab importer',''}, ...
  };
  for index=1:length(formats) % the default loaders are addded after the INI file
    format = formats{index};
    if isempty(format), break; end
    if length(format) < 4, continue; end
    if length(format) < 5, format{5}=''; end
    if length(format) < 6, format{6}=''; end
    % check if format already exists in list
    skip_format=0;
    for j=1:length(loaders)
      this=loaders{j};
      if strcmp(format{1}, this.name)
        skip_format=1;
        break;
      end
    end
    if ~skip_format
      loader.method     = format{1};
      loader.extension  = format{2};
      loader.name       = format{3};
      loader.options    = format{4};
      loader.patterns   = format{5};
      loader.postprocess= format{6};
      loaders = { loaders{:} , loader };
    end
    
  end
  
  config.loaders = loaders; % updated list of loaders

  

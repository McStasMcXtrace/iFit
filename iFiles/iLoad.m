function [data, format] = iLoad(filename, loader)
% [data, loader] = iLoad(file, loader)
%
% imports any data into Matlab
%
% input arguments:
%   file:   file name, or cell of file names, or any Matlab variable, or empty
%             (then popup a file selector)
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
% See also: importdata, load
%
% Part of: iFiles utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. June, 2007.

% calls:    urlread
% optional: uigetfiles, looktxt, unzip, untar, gunzip (can do without)
% private:  iLoad_loader_auto

data = []; format = [];
if nargin == 0, filename=''; end
if nargin < 2,  loader = ''; end

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
  % handle single file name (possibibly with wildcard)
  if ~isempty(find(filename == '*')) | ~isempty(find(filename == '?'))  % wildchar !!#
    [filepath,name,ext,versn]=fileparts(filename);  % 'file' to search
    if isempty(filepath), filepath = pwd; end
    this_dir = dir(filename);
    if isempty(this_dir), this_dir = dir(filepath); end
    index = find(real([this_dir.isdir]) == 0);
    this_dir = char(this_dir.name);
    this_dir = (this_dir(index,:));
    rdir = cellstr(this_dir); % original directory listing as cell
    rdir = strcat([ filepath filesep ], char(rdir{index}));
    filename = cellstr(rdir);
    [data, format] = iLoad(filename, loader);
    return
  end

  [pathstr, name, ext] = fileparts(filename);
  if     strcmp(ext, 'zip'), cmd = 'unzip';
  elseif strcmp(ext, 'tar'), cmd = 'untar';
  elseif strcmp(ext, 'gz'),  cmd = 'gunzip';
  else                       cmd=''; end
  if ~isempty(cmd)
    % this is a compressed file/url. Extract to temporary dir.
    if exist(cmd)
      % extract to temporary dir
      filenames = feval(cmd, filename, tempdir);
      [data, format] = iLoad(filenames, loader); % is now local
      for index=1:length(filenames)
        delete(filenames{index});
      end
      return
    end
  end
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
    [data, format] = iLoad_import(filename, loader);
  end
elseif isempty(filename)
  if exist('uigetfiles')
    [filename, pathname] = uigetfiles;
  else
    [filename, pathname] = uigetfile;
  end
  if isempty(filename), return; end
  filename = strcat(pathname, filesep, filename);
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
      fprintf(1, 'Failed to import file %s with method %s. Ignoring.\n', filename, this_loader.name);
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

  % handle single char loaders
  if ischar(loader)
    tmp=loader; clear loader;
    loader.method = tmp; loader.options='';
  end
  if isempty(loader.method), return; end
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
  return
  
% -----------------------------------------------------------

% private function to make the data pretty looking
function data = iLoad_loader_check(file, data, loader)
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
    data = new_data;
  end
  
  if ~isfield(data, 'Source'),  data.Source = file;         end
  if ~isfield(data, 'Title'),   
    [pathname, filename, ext] = fileparts(file);
    if ~strcmp(loader, 'variable'), data.Title  = [ 'File ' filename ext ' ' name  ];
    else data.Title  = [ 'File ' filename ext ]; end
  end
  if ~isfield(data, 'Date'),    data.Date   = datestr(now); end
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
function loaders = iLoad_loader_auto(file)
  loaders      = {};
  % read user list of loaders which are put 
  if exist('iData_load_ini')
    user_loaders = iData_load_ini;
    loaders = { user_loaders{:} , loaders{:} }; % returns a cell of format specifications (structures)
  end
  
  % read default loaders: method, ext, name, options
  % default importers, when no user specification is given. 
  % These do not have any pattern recognition or postprocess
  formats = {...
    { 'looktxt', '',    'Data (text format with fastest import method)',    '--headers --binary --fast --comment= '}, ...
    { 'looktxt', '',    'Data (text format with fast import method)',       '--headers --binary --comment= '}, ...
    { 'looktxt', '',    'Data (text format)',                               '--headers --comment= '}, ...
    { 'importdata','',  'Matlab importer',''}, ...
    { 'wk1read', 'wk1', 'Lotus1-2-3 (first spreadsheet)',''}, ...
    { 'auread',  'au',  'NeXT/SUN (.au) sound',''}, ...
    { 'wavread', 'wav'  'Microsoft WAVE (.wav) sound',''}, ...
    { 'aviread', 'avi', 'Audio/Video Interleaved (AVI) ',''}, ...
    { 'cdfread', {'nc','cdf'}, 'NetCDF',''}, ...
    { 'netcdf',  {'nc','cdf'}, 'NetCDF 1.0',''}, ...
    { 'fitsread','fits','FITS',''}, ...
    { 'xlsread', 'xls', 'Microsoft Excel (first spreadsheet)',''}, ...
    { 'imread',  {'bmp','jpg','jpeg','tiff','png','ico'}, 'Image',''}, ...
    { 'hdfread', 'h4',  'HDF4',''}, ...
    { 'hdf5read',{'hdf','h5'}, 'HDF5',''}, ...
    { 'load',    'mat', 'Matlab workspace',''}, ...
    { 'csvread', 'csv', 'Comma Separated Values',''}, ...
    { 'dlmread', 'dlm', 'Numerical single block',''}, ...
    { 'xmlread', 'xml', 'XML',''}, ...
  };
  for index=1:length(formats) % the default loaders are addded after the INI file
    format = formats{index};
    if isempty(format), break; end
    loader.method     = format{1};
    loader.extension  = format{2};
    loader.name       = format{3};
    loader.options    = format{4};
    loader.patterns   = '';
    loader.postprocess= '';
    loaders = { loaders{:} , loader };
  end
    
  % read start of file
  fid = fopen(file, 'r');
  if fid == -1
    error([ 'Could not open file ' file ' for reading. Check existence/permissions.' ]);
  end
  file_start = fread(fid, 100000, 'uint8=>char')';
  fclose(fid);
  % loop to test each format for patterns
  formats = loaders;
  loaders={};
  loaders_count=0;
  % identify by extensions
  [duummy, dummy, fext] = fileparts(file);
  for index=1:length(formats)
    loader = formats{index};
    if ~isfield(loader,'extension'), ext=''; else
    ext=loader.extension; end
    if ischar(ext), ext=cellstr(ext); end
    for index_ext=1:length(ext(:))
      if ~isempty(strfind(fext, ext{index_ext}))
        loaders_count = loaders_count+1;
        loaders{loaders_count} = loader;
        break;
      end
    end
  end
  % identify by patterns
  for index=1:length(formats)
    loader = formats{index};
    if ~isstruct(loader), break; end

    if exist(loader.method)
      if length(find(file_start >= 32 & file_start < 127))/length(file_start) < 0.9 ... 
        & strcmp(loader.method, 'looktxt') 
        patterns_found  = 0;
        continue;  % does not use looktxt for binary data files
      end
      patterns_found  = 1;
      if isempty(loader.patterns)  % no pattern to search, just try loader
        patterns_found  = 1;
      else  % check patterns
        if ischar(loader.patterns), loader.patterns=cellstr(loader.patterns); end
        for index_pat=1:length(loader.patterns(:))
          if isempty(strfind(file_start, loader.patterns{index_pat}))
            patterns_found=0;
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
  

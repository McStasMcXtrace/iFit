function out = load(a, varargin)
% d = load(s, file, loader): iData file loader
%
%   @iData/load: imports any data into Matlab/iData object(s)
%   The input argument 'file' should be a file name, or a cell of file names, 
%     or any Matlab variable, or empty (then popup a file selector).
%   The choice of the data file importer is set by default to automatic, so
%     that most common data importers are tested until one works. User may configure
%     a list of prefered loader definitions in a file called iData_load_ini.
%   The optional 3rd argument can be force to use a specific loader list (see below)
%   The input iData object is updated if no output argument is specified.
%
%   Default supported formats include: any text based including CSV, Lotus1-2-3, SUN sound, 
%     WAV sound, AVI movie, NetCDF, FITS, XLS, BMP GIF JPEG TIFF PNG ICO images,
%     HDF4, HDF5, MAT workspace, XML
%   Other specialized formats include: McStas, ILL, SPEC, ISIS/SPE, INX, EDF
%   Compressed files are also supported, with on-the-fly extraction (zip, gz, tar, Z).
%   Distant files are supported through e.g. URLs such as 
%     file://, ftp://, http:// and https://
%   File names may end with an internal anchor reference '#anchor", as used in HTML 
%     links, in which case the members matching the anchor are returned.
%
% input:  s: object or array (iData)
%         file: file name(s) to import (char/cellstr)
%         loader: optional loader method specification (char/struct/cellstr/array of struct)
%               loader = 'auto' (default) test all known data readers until one works
%               loader = 'gui'  manually ask user for the loader(s) to use
%             OR a function name to use as import routine, OR a struct/cell of struct with:
%               loader.method     = function to read data file (char/function_handle)
%               loader.options    = options (char/cell which is then expanded)
%               loader.postprocess= function to act on freshly imported iData (char/function_handle)
%
% The loading process calls first
%   data =loader.method(filename, options...)
% then build the iData object, and optionally calls
%   iData=loader.postprocess(iData)
%
% output: d: single object or array (iData)
% ex:     load(iData,'file'); load(iData); load(iData, 'file', 'gui'); load(a,'','looktxt')
%         load(iData, 'http://file.gz#Data')
%
% Version: $Revision: 1.21 $
% See also: iLoad, save, iData/saveas, iData_load_ini

% calls private/iLoad
% iLoad returns nearly an iData structure.
% inline: load_check_struct, load_clean_metadata
% EF 23/09/07 iData implementation

[files, loaders] = iLoad(varargin{:}); 

if isempty(files), out=[]; return; end
if ~iscell(files),   files = { files }; end
if isstruct(loaders), loaders = { loaders }; end
out = [];
for i=1:length(files)
  filename = '';
  if length(varargin) >= 1 && ischar(varargin{1}), filename = varargin{1}; end
  if isempty(filename) && isstruct(files{i}) && isfield(files{i},'Filename'), filename = files{i}.Filename; end
  if isempty(filename) && isstruct(files{i}) && isfield(files{i},'filename'), filename = files{i}.filename; end
  if isempty(filename) && isstruct(files{i}) && isfield(files{i},'Source'), filename = files{i}.Source; end
  files{i} = load_check_struct(files{i}, loaders, filename);
  this_iData =  iData(files{i});	% convert file content from iLoad into iData
  % specific adjustments for looktxt (default import method)
  [pathname,filename,ext] = fileparts(files{i}.Source);
  try % create MetaData alias if present in structure
    c = this_iData.Data.MetaData; clear c;
    this_iData=setalias(this_iData, 'MetaData', 'Data.MetaData', [ 'MetaData from ' filename ext ]);
    this_iData=load_clean_metadata(this_iData);
  end

  if isfield(files{i},'Headers')
    this_iData.Data.Headers = files{i}.Headers;
    this_iData=setalias(this_iData, 'Headers', 'Data.Headers', [ 'Headers from ' filename ext ]);
  end
  if ~isfield(loaders{i}, 'postprocess')
    loaders{i}.postprocess='';
  end
  if ~isempty(loaders{i}.postprocess)
    % removes warnings
    try
      warn.set = warning('off','iData:setaxis');
      warn.get = warning('off','iData:getaxis');
      warn.get = warning('off','iData:get');
    catch
      warn = warning('off');
    end
    if ~iscellstr(loaders{i}.postprocess)
      loaders{i}.postprocess = cellstr(loaders{i}.postprocess);
    end
    % apply post-load routine: this may generate more data sets
    for j=1:length(loaders{i}.postprocess)
      if ~isempty(loaders{i}.postprocess{j})
        this_iData = feval(loaders{i}.postprocess{j}, this_iData);
      end
    end
    % reset warnings
    try
      warning(warn.set);
      warning(warn.get);
    catch
      warning(warn);
    end
  elseif ~isempty(loaders{i}.postprocess)
    iData_private_warning(mfilename,['Can not find post-process function ' loaders{i}.postprocess ' for data format ' loaders{i}.name ]);
  end
  out = [ out this_iData ];
end %for i=1:length(files)
for i=1:length(out)
  out(i).Command={[ out(i).Tag '=load(iData,''' out(i).Source ''');' ]};
  %this_iData = iData_private_history(this_iData, mfilename, a, files{i}.Source);
end % for

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),out);
end

% ------------------------------------------------------------------------------
function a=load_clean_metadata(a, loaders, filenames)
% test each field of MetaData and search for equal aliases
  this = a.Data.MetaData;
  meta_names = fieldnames(this);
  alias_names=getalias(a);
  %treat each MetaData
  for index=1:length(meta_names)
    if numel(getfield(this, meta_names{index})) > 1000
    for index_alias=1:length(alias_names)
      % is it big and equal to an alias value ?
      if isequal(getfield(this, meta_names{index}), get(a, alias_names{index_alias}))
        % yes: store alias in place of MetaData
        this = setfield(this, meta_names{index}, getalias(a, alias_names{index_alias}));
        break
      end
    end % for
    end % if
  end
  a.Data.MetaData = this;

% ------------------------------------------------------------------------------
function s=load_check_struct(data, loaders, filename)
  if nargin < 3, filename=''; end
  if isempty(filename), filename=pwd; end
  if iscell(filename),  filename=filename{1}; end
  
  % transfer some standard fields as possible
  if ~isstruct(data)          s.Data = data; else s=data; end
  if isfield(data, 'Source'), s.Source = data.Source; 
  else                        s.Source = filename; end
  if isfield(data, 'Title'),  s.Title = data.Title; 
  else 
    [pathname, filename, ext] = fileparts(filename);
    s.Title  = [ 'File ' filename ext ];
  end
  if isfield(data, 'Date'),   s.Date = data.Date; 
  else                        s.Date   = datestr(now); end
  if isfield(data, 'Label'),  s.Label = data.Label; end
  if ~isfield(s, 'Format'),
    s.Format  = loaders{1}.name; 
  end

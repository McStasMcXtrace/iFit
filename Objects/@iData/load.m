function out = load(a, varargin)
% d = load(s, file, loader, ...): iData file loader
%
%   @iData/load: imports any data into Matlab/iData object(s)
%   The input argument 'file' should be a file name, or a cell of file names, 
%     or any Matlab variable, or empty (then popup a file selector).
%   The choice of the data file importer is set by default to automatic, so
%     that most common data importers are tested until one works. User may configure
%     a list of prefered loader definitions in a file called iData_load_ini.
%   The optional 3rd argument can be set to use a specific loader list (see below)
%     load(iData, filename, loader)
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
%  Type <a href="matlab:doc(iData,'Load')">doc(iData,'Load')</a> to access the iFit/Load Documentation.
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
%         additional arguments are passed to the import routine.
%
% The loading process calls first
%   data =loader.method(filename, options...)
% then build the iData object, and optionally calls
%   iData=loader.postprocess(iData)
%
% output: d: single object or array (iData)
% ex:     load(iData,'file'); load(iData); load(iData, 'file', 'gui'); load(a,'','looktxt')
%         load(iData, [ ifitpath 'Data/peaks.hdf5' ], 'HDF')
%         load(iData, 'http://file.gz#Data')
%
% Version: $Revision: 1.40 $
% See also: iLoad, save, iData/saveas, Loaders

% calls private/iLoad
% iLoad returns nearly an iData structure.
% inline: load_check_struct, load_clean_metadata
% EF 23/09/07 iData implementation

[files, loaders] = iLoad(varargin{:}); % import files as structures HERE

if isempty(files), out=[]; return; end
if isstruct(files) && length(files) == 1 && isfield(files,'loaders')
    out=files;
    return;
end
if ~iscell(files),   files   = { files }; end
if ~iscell(loaders), loaders = { loaders }; end
out = [];

for i=1:numel(files)
  filename = '';
  if length(varargin) >= 1 && ischar(varargin{1}), filename = varargin{1}; end
  if isempty(filename) && isstruct(files{i}) && isfield(files{i},'Filename'), filename = files{i}.Filename; end
  if isempty(filename) && isstruct(files{i}) && isfield(files{i},'filename'), filename = files{i}.filename; end
  if isempty(filename) && isstruct(files{i}) && isfield(files{i},'Source'), filename = files{i}.Source; end
  files{i} = load_check_struct(files{i}, loaders, filename);
  if isfield(files{i},'Data') && isstruct(files{i}.Data) && any(cellfun('isclass', struct2cell(files{i}.Data), 'iData'))
    this_iData = [];
    struct_data = struct2cell(files{i}.Data);
    files{i} = {};  % free memory
    for index=1:length(struct_data)
      if isa(struct_data{index}, 'iData')
        this_iData = [ this_iData struct_data{index} ];
      end
    end
  else
    this_iData =  iData(files{i});	% convert file content from iLoad into iData
    if ~isempty(this_iData)
      % specific adjustments for looktxt (default import method)
      [pathname,filename,ext] = fileparts(files{i}.Source);
      if isfield(this_iData.Data, 'MetaData')
        this_iData=setalias(this_iData, 'MetaData', 'Data.MetaData', [ 'MetaData from ' filename ext ]);
        this_iData=load_clean_metadata(this_iData);
      end
      
      if isempty(loaders) || ~isfield(loaders{i}, 'postprocess')
        loaders{i}.postprocess='';
      end
      if isempty(loaders{i}.postprocess) && isfield(files{i},'Loader')
        if isfield(files{i}.Loader, 'postprocess')
          loaders{i}.postprocess = files{i}.Loader.postprocess;
        end
      end
      files{i} = {};  % free memory
      if ~isempty(loaders{i}.postprocess)
        % removes warnings
        iData_private_warning('enter',mfilename);
        if ~iscellstr(loaders{i}.postprocess)
          loaders{i}.postprocess = cellstr(loaders{i}.postprocess);
        end
        % apply post-load routine: this may generate more data sets
        for j=1:length(loaders{i}.postprocess)
          if ~isempty(loaders{i}.postprocess{j})
            % disp([ mfilename ': Calling post-process ' loaders{i}.postprocess{j} ]);
            try
              this_iData = feval(loaders{i}.postprocess{j}, this_iData);
              this_iData = setalias(this_iData, 'postprocess', loaders{i}.postprocess{j});
              this_iData = iData_private_history(this_iData, loaders{i}.postprocess{j}, this_iData);
            catch
              iData_private_warning(mfilename, [ 'Error when calling post-process ' loaders{i}.postprocess{j} ]);
            end
          end
        end
        % reset warnings
        iData_private_warning('exit',mfilename);
      elseif ~isempty(loaders{i}.postprocess)
        iData_private_warning(mfilename,['Can not find post-process function ' loaders{i}.postprocess ' for data format ' loaders{i}.name ]);
      end
    end
  end
  out = [ out this_iData ];
end %for i=1:length(files)

% clean 'out' for unique entries (e.g. mcstas/mcstas.sim creates duplicates)
[b, i1] = unique(get(out, 'Source')); % name is unique
if 1 < length(out) && length(i1) < length(out)
  % some data sets seem to be duplicated: make additional tests
  [b, i2] = unique(sum(out, 0));      % sum of Signal is unique
  [b, i3] = unique(get(out, 'Title'));% Title is unique
  [b, i4] = unique(get(out, 'Label'));% Title is unique
  i5 = get(out,'Signal');
  if iscell(i5)
    [b, i5] = unique(cellfun('prodofsize',i5)); % size of Signal is unique
  else
    i5=[];
  end
  i = unique([i1(:) ; i2(:) ; i3(:) ; i4(:) ; i5(:) ]);
  out = out(i);
end

for i=1:numel(out)
  out(i).Command{end+1}=[ out(i).Tag '=load(iData,''' out(i).Source ''');' ];
  if isempty(out(i).DisplayName)
    [p,f,e] = fileparts(filename);
    out(i) = set(out(i),'DisplayName',[ f e ]);
  end
  %this_iData = iData_private_history(this_iData, mfilename, a, files{i}.Source);
end % for

if nargout == 0 && ~isempty(inputname(1))
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
  if isfield(data, 'Date'),   s.Date   = data.Date; 
  else                        s.Date   = clock; end
  if isfield(data, 'Label'),  s.Label = data.Label; end
  if ~isfield(s, 'Format'),
    s.Format  = loaders{1}.name; 
  end
  

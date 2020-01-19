function out = load(a, varargin)
% LOAD Load data from file(s) into workspace.
%   s = LOAD(estruct, file) loads the given file as an object.
%   The input argument 'file' should be a file name, or a cell of file names, 
%   or any Matlab variable, or empty (then popup a file selector).
%   The choice of the data file importer is set by default to automatic, so
%   that most common data importers are tested until one works. 
%   The list of preferred loader definitions can be set into an iLoad_ini file or
%   iLoad.ini in the preferences directory. Refer to 'iload' for more information.
%
%   s = LOAD(estruct, file, loader) same as above, and specifies the loader to
%   use. The loader can be an exact function name, such as 'read_hdf5'
%   or a keyword to search into all loaders, such as e.g. 'hdf5'.
%   OR a struct/cell of struct with:
%               loader.method     = function to read data file (char/function_handle)
%               loader.options    = options (char/cell which is then expanded)
%               loader.postprocess= char/function_handle to act on freshly imported estruct
%   The loading process calls first
%     data =loader.method(filename, options...)
%   then build the estruct object, and optionally calls
%     estruct=loader.postprocess(estruct)
%
%   s = LOAD(estruct, file, loader, ...) pass additional arguments to the loader.
%
%   LOAD(estruct, file, 'auto') is equivalent to LOAD(estruct, file).
%
%   LOAD(estruct, file, 'gui') displays a list dialogue to select a loader to use.
%
%   LOAD(estruct, ...) loads data, and updates initial estruct object.
%   The input estruct object is updated if no output argument is specified.
%
%   Default supported formats include: any text based including CSV, Lotus1-2-3, SUN sound, 
%     WAV sound, AVI movie, NetCDF, FITS, XLS, BMP GIF JPEG TIFF PNG ICO images,
%     HDF4, HDF5, MAT workspace, XML, CDF, JSON, YAML, IDL
%   Other specialized formats include: McStas, ILL, SPEC, ISIS/SPE, INX, EDF, Mantid.
%     SIF, MCCD/TIFF, ADSC, CBF, Analyze, NifTI, STL,PLY,OFF, CIF/CFL,
%     EZD/CCP4, Bruker Varian and JEOL NMR formats, Bruker OPUS, LabView LVM and TDMS
%     Agilent and Thermo Finnigan MS, Quantum Design VMS, ENDF
%   Compressed files are also supported, with on-the-fly extraction (zip, gz, tar, Z).
%   Distant files are supported through e.g. URLs such as 
%     file://, ftp://, http:// and https://
%   File names may end with an internal anchor reference '#anchor", as used in HTML 
%     links, in which case the members matching the anchor are returned.
%   The list of supported formats obtained with:
%     load(estruct,'formats') or iLoad('formats')
%
%  Type <a href="matlab:doc(estruct,'Load')">doc(estruct,'Load')</a> to access the iFit/Load Documentation.

% Example: isa('estruct', load(estruct, 'peaks.hdf5'))
% Version: $Date$ $Version$ $Author$
% See also: iLoad, save, estruct/saveas, Loaders

[files, loaders] = iLoad(varargin{:}); % import files as structures HERE

if isempty(files), out=[]; return; end
if isstruct(files) && length(files) == 1 && isfield(files,'loaders')
    out=files;
    return;
end
% convert struct array to cell array, if true
if isstruct(files) && numel(files) > 1 && numel(loaders) == 1
  new_files = cell(1,numel(files));
  new_loaders = new_files;
  for index=1:numel(files)
    new_files{index} = files(index);
    new_loaders{index} = loaders;
  end
  files = new_files; new_files = [];
  loaders = new_loaders; new_loaders = [];
end
if ~iscell(files),   files   = { files }; end
if ~iscell(loaders), loaders = { loaders }; end
out = [];
loader = [];
for i=1:numel(files)
  filename = '';
  if isempty(files{i}), continue; end
  if length(varargin) >= 1 && ischar(varargin{1}), filename = varargin{1}; end
  if i <= numel(loaders), loader = loaders{i}; end
  
  this_estruct = load_single_file(files{i}, loader, filename);
  
  % transpose to make columns
  if numel(out) > 1 && size(out,2) > 1 && size(out,1) == 1, out=out'; end
  if numel(this_estruct) > 1 && size(this_estruct,2) > 1 && size(this_estruct,1) == 1, this_estruct=this_estruct'; end
  out = [ out ; this_estruct ];
  clear this_estruct
end %for i=1:length(files)

% clean 'out' for unique entries (e.g. mcstas/mcstas.sim creates duplicates)
[b, i1] = unique(get(out, 'Source')); % name is unique
if 1 < length(out) && length(i1) < length(out)
  % some data sets seem to be duplicated: make additional tests
  % look for similarities
  sources = get(out, 'Source');
  titls   = get(out, 'Title');
  labs    = get(out, 'Label');
  sums    = sum(out, 0); % total signal
  i       = 1:length(out);
  for index=1:length(out)
    j = find(strcmp(sources{index}, sources(:)) & strcmp(titls{index}, titls(:)) & strcmp(labs{index}, labs(:)) & sums(index) == sums(:));
    if length(j) > 1, i(j(2:end)) = 0; end
  end
  removed = find(i == 0);
  i       = unique(i(i>0));
  if length(out) > length(i)
    warning('%s: Removing duplicated data sets %i -> %i', mfilename, length(out), length(i))
    warning([ mfilename ': ' char(out(removed)) ])
    out = out(i);
  end
end

for i=1:numel(out)
  out(i).Command{end+1}=[ out(i).Tag '=load(estruct,''' out(i).Source ''');' ];
  if isempty(out(i).DisplayName) && isempty(out(i).Label)
    [p,f,e] = fileparts(out(i).Source);
    out(i) = set(out(i),'Label',[ f e ]);
  end
  %this_estruct = estruct_private_history(this_estruct, mfilename, a, files{i}.Source);
end % for

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),out);
end

% ----------------------------------------------------------------------
function s=load_check_struct(data, loaders, filename)
% check final structure and add missing fields
  if nargin < 3, filename=''; end
  if isempty(filename), filename=pwd; end
  if iscell(filename),  filename=filename{1}; end
  
  if isstruct(data) && numel(data) > 1
      s = [];
      for index=1:numel(data)
          s = [ s ; load_check_struct(data(index), loaders, filename) ];
      end
      return
  end

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
    if iscell(loaders) s.Format  = loaders{1}.name; 
    else               s.Format  = loaders.name; 
    end
  end

% ----------------------------------------------------------------------
function this = load_eval_postprocess(this, postprocess)
% evaluate the postprocess in a reduced environment, with only 'this'
  this0 = this;
  try
    % disp([ mfilename ': Calling post-process '  postprocess ])
    if isvarname(postprocess) && exist(postprocess) == 2
      this = feval(postprocess, this);
    elseif isempty(postprocess == '=')
      this = eval(postprocess);
    else
      eval(postprocess);
    end
    
    if ~isa(this, 'estruct'), this = this0;
    else
      this = setalias(this, 'postprocess', postprocess);
      this = history(this, postprocess, this);
    end
  catch ME
    warning(getReport(ME));
    warning([mfilename ': Error when calling post-process ' postprocess '. file: ' this.Source ]);
  end
  
% ------------------------------------------------------------------------------
function this_estruct = load_single_file(file, loader, filename)
    
  this_estruct = [];
  
  % handle array of struct
  if numel(file) > 1
    for index=1:numel(file)
      if isstruct(file)
        this_estruct = [ this_estruct ; load_single_file(file(index), loader, filename) ];
      elseif iscell(file)
        this_estruct = [ this_estruct ; load_single_file(file{index}, loader, filename) ];
      end
    end
    return
  end
  
  if isstruct(file) && isempty(filename)
    f = fieldnames(file);
    index = [ find(strcmpi(f,'filename'),1) ;find(strcmpi(f,'file_name'),1) ;find(strcmpi(f,'source'),1) ];
    if ~isempty(index)
      filename = file.(f{index(1)}); 
    end
  end
  
  % check the returned iLoad structure
  file = load_check_struct(file, loader, filename);
  if isfield(file,'Data') && isstruct(file.Data) && any(cellfun('isclass', struct2cell(file.Data), 'estruct'))
    % a structure containing an estruct
    this_estruct = [];
    struct_data = struct2cell(file.Data);
    file = {};  % free memory
    for index=1:length(struct_data)
      if isa(struct_data{index}, 'estruct')
        if isa(struct_data{index}.Data, 'uint8')
          struct_data{index}.Data = hlp_deserialize(struct_data{index}.Data);
        end
        this_estruct = [ this_estruct struct_data{index} ];
        struct_data{index} = '';
      end
    end
    clear struct_data
  else
    % usually a structure from iLoad

    % convert file content from iLoad into estruct
    this_estruct = struct2estruct(file);	
    % assign default Signal and axes
    this_estruct = axescheck(this_estruct);
    % post-processing
    if ~isempty(this_estruct)
      if isempty(loader) || ~isfield(loader, 'postprocess')
        if ischar(loader), name=loader; loader=[]; loader.name=name; end
        loader.postprocess='';
      end
      if isempty(loader.postprocess) && isfield(file,'Loader')
        if isfield(file.Loader, 'postprocess')
          loader.postprocess = file.Loader.postprocess;
        end
      end
      file = {};  % free memory
      if ~isempty(loader.postprocess)
        if ~iscell(loader.postprocess)
          loader.postprocess = cellstr(loader.postprocess);
        end
        % apply post-load routine: this may generate more data sets
        for j=1:length(loader.postprocess)
          if ~isempty(loader.postprocess{j})
            % call private method (see below)
            this_estruct = load_eval_postprocess(this_estruct, loader.postprocess{j});
          end
        end
        % reset warnings
      elseif ~isempty(loader.postprocess)
        warning([ mfilename ': Can not find post-process function ' loader.postprocess ' for data format ' loaders{i}.name ]);
      end
    end
  end

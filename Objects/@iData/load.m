function out = load(a, varargin)
% LOAD Load data from file(s) into workspace.
%   s = LOAD(iData, file) loads the given file as an object.
%   The input argument 'file' should be a file name, or a cell of file names, 
%   or any Matlab variable, or empty (then popup a file selector).
%   The choice of the data file importer is set by default to automatic, so
%   that most common data importers are tested until one works. 
%   The list of preferred loader definitions can be set into an iLoad_ini file or
%   iLoad.ini in the preferences directory. Refer to 'iload' for more information.
%   As opposed to IDATA(file), the LOAD method can specify a given loader, and
%   applies further post-processing filters.
%
%   s = LOAD(iData, file, loader) same as above, and specifies the loader to
%   use. The loader can be an exact function name, such as 'read_hdf5'
%   or a keyword to search into all loaders, such as e.g. 'hdf5'.
%   OR a struct/cell of struct with:
%               loader.method     = function to read data file (char/function_handle)
%               loader.options    = options (char/cell which is then expanded)
%               loader.postprocess= char/function_handle to act on freshly imported iData
%   The loading process calls first
%     data =loader.method(filename, options...)
%   then build the iData object, and optionally calls
%     iData=loader.postprocess(iData)
%
%   s = LOAD(iData, file, loader, ...) pass additional arguments to the loader.
%
%   LOAD(iData, file, 'auto') is equivalent to LOAD(iData, file).
%
%   LOAD(iData, file, 'gui') displays a list dialogue to select a loader to use.
%
%   LOAD(iData, ...) loads data, and updates initial iData object.
%   The data can be anything, including a structure, cell, array, ...
%   The input iData object is updated if no output argument is specified.
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
%     load(iData,'formats') or iLoad('formats')
%
%  Type <a href="matlab:doc(iData,'Load')">doc(iData,'Load')</a> to access the iFit/Load Documentation.
%
% Example: isa(load(iData, 'peaks.hdf5'),'iData')
% Version: $Date$ $Version$ $Author$
% See also: iLoad, save, iData/saveas, Loaders

% check for 'no postprocess/raw' options
use_post_process = true;
varg = varargin;
for index=1:numel(varargin)
  if ischar(varargin{index}) && any(strcmp(varargin{index},{'raw','no postprocess'}))
    use_post_process = false;
    varg{index} = [];
  end
end
varargin = varg;

[files, loaders] = iLoad(varargin{:}); % import files as structures HERE

if isempty(files), out=[]; return; end
if isstruct(files) && length(files) == 1 && isfield(files,'loaders')
    out=files;
    return;
end

% convert struct array to cell array, if true
if isstruct(files) && numel(files) > 1 && numel(loaders) == 1
  new_files   = cell(1,numel(files));
  new_loaders = new_files;
  for index=1:numel(files)
    new_files{index}   = files(index);
    new_loaders{index} = loaders;
  end
  files   = new_files;   new_files = [];
  loaders = new_loaders; new_loaders = [];
end
if ~iscell(files),   files   = { files }; end
if ~iscell(loaders), loaders = { loaders }; end
out    = [];
loader = [];

% convert iLoad structs to iData (no post-process yet)
for i=1:numel(files)
  if isempty(files{i}), continue; end

  this_iData = struct2iData(files{i}); % convert into new iData
  
  % and call any postprocess (if any)
  if use_post_process && isfield(this_iData, 'postprocess')
    if a.verbose
      try; disp([ mfilename ': calling post-process ' char(this_iData.postprocess) ' for ' char(this_iData,'short') ]); end
    end
    this_iData = private_postprocess(this_iData, this_iData.postprocess); % can return an array
    if numel(this_iData) > 1
      % need to check for duplicates (when post-process creates new data sets)
      try
        this_iData = private_remove_duplicates(this_iData);
      end
    end
  end

  out = [ out ; this_iData(:) ];
  clear this_iData
end %for i=1:length(files)

% set Command history and Label/DisplayName
history(out, mfilename, a, varargin{:});
for i=1:numel(out)
  if isempty(out(i).DisplayName) && isempty(out(i).Label)
    [p,f,e] = fileparts(out(i).Source);
    out(i)  = set(out(i),'Label',[ f e ]);
  end
end % for

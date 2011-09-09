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
%     file://, ftp:// and http://
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
%
% Version: $Revision: 1.19 $
% See also: iLoad, save, iData/saveas, iData_load_ini

% calls private/iLoad
% iLoad returns nearly an iData structure.
% EF 23/09/07 iData implementation

if isempty(varargin), [files, loaders] = iLoad; 
else [files, loaders] = iLoad(varargin{:}); end

if isstruct(files),   files = { files }; end
if isstruct(loaders), loaders = { loaders }; end
out = [];
for i=1:length(files)
  this_iData =  iData(files{i});	% convert structure from iLoad into iData
  % specific adjustments for looktxt (default import method)
  [pathname,filename,ext] = fileparts(files{i}.Source);
  try % create MetaData alias if present in structure
    c = this_iData.Data.MetaData; clear c;
    this_iData=setalias(this_iData, 'MetaData', 'Data.MetaData', [ 'MetaData from ' filename ext ]);
  catch
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


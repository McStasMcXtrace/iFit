function out = load(a, varargin)
% d=load(s, file, loader): iData file loader
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
  this_iData =  iData(files{i});
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
  
  if ~isempty(loaders{i}.postprocess)
    this_iData = feval(loaders{i}.postprocess, this_iData);
  end
  
  this_iData = iData_private_history(this_iData, mfilename, a, files{i}.Source);
  out = [ out this_iData ];  
end

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),out);
end


function out = load(a, varargin)
% d=load(s, file, ...): iData file loader
%
%   @iData/load: imports any data into Matlab/iData object(s)
%   The input arguments 'file' should be a file name, or a cell of file names, 
%   or any Matlab variable, or empty (then popup a file selector).
%   The input iData object is updated if no output argument is specified.
%
% input:  s: object or array (iData)
%         file: file name to import (char)
% output: d: single object or array (iData)
% ex:     load(iData,'file') or load(iData)
%
% See also: iLoad, save, iData/saveas

% calls iFiles/iLoad
% iLoad returns nearly an iData structure.
% EF 23/09/07 iData implementation

if isempty(varargin), files = iLoad; 
else files = iLoad(varargin); end

files = files{1};
if isstruct(files), files = { files }; end
out = [];
for i=1:length(files)
  out = [ out iData(files{i}) ];
  % specific adjustments for looktxt (default import method)
  [pathname,filename,ext] = fileparts(files{i}.Source);
  try % create MetaData alias if present in structure
    c = out(i).Data.MetaData; clear c;
    out(i)=setalias(out(i), 'MetaData', 'Data.MetaData', [ 'MetaData from ' filename ext ]);
  catch
  end

  if isfield(files{i},'Headers')
    out(i).Data.Headers = files{i}.Headers;
    out(i)=setalias(out(i), 'Headers', 'Data.Headers', [ 'Headers from ' filename ext ]);
  end
  
  out(i) = iData_private_history(out(i), mfilename, a, files{i}.Source);
  
end

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),out);
end


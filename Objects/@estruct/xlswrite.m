function f=xlswrite(a, filename, options)
% XLSWRITE  Write object to Microsoft Excel spreadsheet file.
%   F = XLSWRITE(S) creates a XLS file named from the object Tag.
%   Requires Microsoft Excel to be installed.
%
%   F = XLSWRITE(S, FILENAME) saves object into given XLS file.
%
% Example: f=xlswrite(estruct(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/saveas

  if nargin < 2, filename = []; end
  if nargin < 3, options = ''; end
 
  f=saveas(a, filename, 'xls',options);

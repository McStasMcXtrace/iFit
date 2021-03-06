function f=csvwrite(a, filename)
% CSVWRITE Write a comma-separated value file.
%   F = CSVWRITE(S, FILENAME) writes object S into FILENAME as comma-separated 
%   values.
%
% Example: f=csvwrite(iData(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/saveas


  if nargin < 2, filename = []; end
 
  f=saveas(a, filename, 'csv');


function f=cdfwrite(a, filename)
% CDFWRITE Write data to a CDF file.
%   F = CDFWRITE(S) creates a CDF file named from the object Tag.
%
%   F = CDFWRITE(S, FILENAME) saves object into given CDF file.
%
% Example: f=cdfwrite(iData(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/saveas

  if nargin < 2, filename = []; end
  if nargin < 3, options = ''; end
 
  f=saveas(a, filename, 'cdf');

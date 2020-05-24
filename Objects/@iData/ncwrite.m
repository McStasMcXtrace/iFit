function f=ncwrite(a, filename, options)
% NCWRITE Write object to NetCDF file.
%   F = NCWRITE(S) creates a NetCDF file named from the object Tag.
%
%   F = NCWRITE(S, FILENAME) saves object into given NetCDF file.
%
% Example: f=ncwrite(iData(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/saveas

  if nargin < 2, filename = []; end
  if nargin < 3, options = ''; end
 
  f=saveas(a, filename, 'nc');

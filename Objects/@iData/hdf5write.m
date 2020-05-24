function f=hdf5write(a, filename)
% HDF5WRITE Write object to HDF5 file.
%   F = HDF5WRITE(S) creates a HDF5 file named from the object Tag.
%
%   F = HDF5WRITE(S, FILENAME) saves object into given HDF5 file.
%
% Example: f=hdf5write(iData(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/saveas

  if nargin < 2, filename = []; end
  if nargin < 3, options = ''; end
 
  f=saveas(a, filename, 'hdf5');

function f=fitswrite(a, filename)
% FITSWRITE Write 2D object to FITS file.
%   F = FITSWRITE(S) creates a FITS file named from the object Tag.
%
%   F = FITSWRITE(S, FILENAME) saves object into given FITS file.
%
% Example: f=fitsfile(estruct(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/saveas

  if nargin < 2, filename = []; end
 
  f=saveas(a, filename, 'fits');

function f=imwrite(a, filename, format)
% IMWRITE Write object to graphics file.
%   IMWRITE(A,FILENAME,FMT) writes the object A to the file specified by
%   FILENAME in the format specified by FMT.
%
%   Supported formats are:
%     GIF BMP PBM PCX PGM PNM PPM RAS XWD HDF4 TIFF PNG ART(ascii)
%     Default is PNG.
%
% Example: f=imwrite(iData(peaks)); ~isempty(dir(f)); delete(f);
% Version: $Date$ $Version$ $Author$
% See also iData, iData/saveas

  if nargin < 2, filename = []; end
  if nargin < 3, format = 'png'; end
  
  f = saveas(a, filename, format);

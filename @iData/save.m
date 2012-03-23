function [filename,format] = save(a, varargin)
% f = save(s, filename, format, options) : save iData object into various data formats
%
%   @iData/save function to save data sets
%     This function saves the content of iData objects. The default format is 'm'.
%   save(iData,'formats')
%     prints a list of supported export formats.
%   save(iData,'file.ext')            determine file format from the extension
%   save(iData,'file','format')       sets file format explicitly
%   save(iData,'file','format clean') sets file format explicitly and remove NaN and Inf.
%
%  Type <a href="matlab:doc(iData,'Save')">doc(iData,'Save')</a> to access the iFit/Save Documentation.
%
% input:  s: object or array (iData)
%         filename: name of file to save to. Extension, if missing, is appended (char)
%                   If the filename already exists, the file is overwritten.
%                   If given as filename='gui', a file selector pops-up
%         format: data format to use (char), or determined from file name extension
%           'm'    save as a flat Matlab .m file (a function which returns an iData object or structure)
%           'mat'  save as a '.mat' binary file (same as 'save', default)
%           'hdf5' save as an HDF5 data set
%           'nc'   save as NetCDF 
%         as well as other lossy formats
%           'hdf4' save as an HDF4 immage
%           'fig'  save as a Matlab figure
%           'edf'  EDF ESRF format for 1D and 2D data sets
%           'gif','bmp' save as an image (no axes, only for 2D data sets)
%           'png','tiff','jpeg','ps','pdf','ill','eps' save as an image (with axes)
%           'xls'  save as an Excel sheet (requires Excel to be installed)
%           'csv'  save as a comma separated value file
%           'svg'  save as Scalable Vector Graphics (SVG) format
%           'wrl'  save as Virtual Reality VRML 2.0 file
%           'dat'  save as Flat text file with comments
%           'fits' save as FITS binary image (only for 2D objects)
%           'hdr'  save as HDR/IMG Analyze MRI volume (3D)
%           'stl'  save as STL stereolithography (geometry), binary
%           'stla' save as STL stereolithography (geometry), ascii
%           'off'  save as Object File Format (geometry), ascii
%           'ply'  save as PLY (geometry), ascii
%           'gui' when filename extension is not specified, a format list pops-up
%         options: specific format options, which are usually plot options
%           default is 'view2 axis tight'
%
% output: f: filename(s) used to save data (char)
% ex:     b=save(a, 'file', 'm');
%         b=save(a, 'file', 'svg', 'axis tight');
%
% Contributed code (Matlab Central): 
%   plot2svg:   Juerg Schwizer, 22-Jan-2006 
%   iData_private_save_hdfnc
%   pmedf_write
%   fitswrite:  R. G. Abraham, Institute of Astronomy, Cambridge University (1999)
%   stlwrite
%
% Version: $Revision: 1.3 $
% See also iData, iData/saveas, iData/load, iData/getframe, save, saveas

[filename,format] = saveas(a, varargin{:});


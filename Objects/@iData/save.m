function [filename,format] = save(a, varargin)
% f = save(s, filename, format, options) : save iData object into various data formats
%
%   @iData/save function to save data sets
%     This function saves the content of iData objects. The default format is 'm'.
%   save(iData,'formats')
%     prints a list of supported export formats.
%   save(s,'file.ext')            determine file format from the extension
%   save(s,'file','format')       set file format explicitly
%   save(s,'file','format clean') set file format explicitly and remove NaN and Inf.
%   save(s,'file','format data')  save only the 'Data' part of the object. 
%   save(s, get(s,'Source'), 'format')  save with the same base filename(s)
%
%     To load back an object from a m-file, type its file name at the prompt.
%     To load back an object from a mat-file, type 'load filename.mat' at the prompt.
%
%  Type <a href="matlab:doc(iData,'Save')">doc(iData,'Save')</a> to access the iFit/Save Documentation.
%
% input:  s: object or array (iData)
%         filename: name of file to save to. Extension, if missing, is appended (char)
%                   If the filename already exists, the file is overwritten.
%                   If given as filename='gui', a file selector pops-up
%                   If the filename is empty, the object Tag is used.
%         format: data format to use (char), or determined from file name extension
%           'cdf'  save as CDF (not recommended)
%           'hdf5' save as an HDF5 data set ('nxs','n5','h5' also work)
%           'lamp' save as LAMP Processed Workspace, i.e. 'nxs lamp data' (HDF5)
%           'm'    save as a flat Matlab .m file (a function which returns an iData object or structure)
%           'mantid' save as Mantid Processed Workspace, i.e. 'nxs mantid data' (HDF5)
%           'mat'  save as a serialized '.mat' binary file (fast 'save', DEFAULT)
%           'nc'   save as NetCDF
%         as well as other lossy formats
%           'art'  save as ASCII art
%           'avi'  save as an AVI movie
%           'csv'  save as a comma separated value file
%           'dae'  save as Collada model
%           'dat'  save as Flat text file with comments
%           'edf'  EDF ESRF format for 1D and 2D data sets
%           'fig'  save as a Matlab figure
%           'fits' save as FITS binary image (only for 2D objects)
%           'gif','bmp','png','tiff','jpeg' save as an image (no axes, only for 2D data sets)
%           'hdf4' save as an HDF4 image
%           'hdr'  save as HDR/IMG Analyze MRI volume (3D/4D)
%           'html' save as Hypertext Markup Language document, appended to any existing document.
%           'inx'  save as an ILL Inelastic Neutron Scattering data (only 2D neutron)
%           'json' save as JSON JavaScript Object Notation, ascii
%           'kml'  save as KML GoogleEarth model
%           'mrc'  save as MRC map file (3/4D)
%           'nii'  save as NifTi Neuroimaging Informatics Technology Initiative (3/4D)
%           'npy'  save as Numpy binary array
%           'ps','pdf','ill','eps' save as an image (with axes)
%           'ppm','pgm','pbm'
%           'off'  save as Object File Format (geometry), ascii
%           'ply'  save as PLY (geometry), ascii
%           'spe'  save as ISIS SPE (Mslice/Horace)
%           'sqw'  save as McStas SQW Isotropic S(q,w)
%           'stl'  save as STL stereolithography (geometry), binary
%           'stla' save as STL stereolithography (geometry), ascii
%           'svg'  save as Scalable Vector Graphics (SVG) format
%           'vtk'  save as VTK ascii (<1e5 elements) or binary (3/4D)
%           'wrl'  save as Virtual Reality VRML 2.0 file
%           'x3d'  save as X3D (geometry) file, ascii
%           'xhtml' save as embedded HTML/X3D file (using Flash plugin for rendering)
%           'xls'  save as an Excel sheet (requires Excel to be installed)
%           'xml'  save as an XML file, ascii
%           'yaml' save as YAML format, ascii
%
%           'gui' when filename extension is not specified, a format list pops-up
%         options: specific format options, which are usually plot options
%           default is 'view2 axis tight'
%
% output: f: filename(s) used to save data (char)
% ex:     b=save(a, 'file', 'm');
%         b=save(a, 'file', 'svg', 'axis tight');
%         b=save(a, 'file', 'hdf data');
%
% Version: $Date$
% See also iData, iData/saveas, iData/load, iData/getframe, save, saveas,
% iData/plot

[filename,format] = saveas(a, varargin{:});


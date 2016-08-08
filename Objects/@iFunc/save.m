function [filename,format] = save(a, varargin)
% f = save(s, filename, format, options) : save iFunc object into various data formats
%
%   @iFunc/save function to save models
%     This function saves the content of iFunc objects. The default format is 'm'.
%   save(iFunc,'formats')
%     prints a list of supported export formats.
%   save(iFunc,'file.ext')            determine file format from the extension
%   save(iFunc,'file','format')       sets file format explicitly
%     To load back a model from an m-file, type its file name at the prompt.
%     To load back a model from an mat-file, type 'load filename.mat' at the prompt.
%
% input:  s: object or array (iFunc)
%         filename: name of file to save to. Extension, if missing, is appended (char)
%                   If the filename already exists, the file is overwritten.
%                   If given as filename='gui', a file selector pops-up
%         format: data format to use (char), or determined from file name extension
%           'json' save as JSON JavaScript Object Notation, ascii
%           'm'    save as a flat Matlab .m file (a function which returns an iFunc object or structure)
%           'mat'  save as a '.mat' binary file (same as 'save', DEFAULT)
%           'yaml' save as YAML format, ascii
%           'xml'  save as XML file, ascii
%         as well as other lossy formats
%           'fig'  save as a Matlab figure
%           'gif','bmp','png','tiff','jpeg','ps','pdf','ill','eps' save as an image
%           'hdf4' save as an HDF4 immage
%           'html' save as Hypertext Markup Language document, appended to any existing document
%
%           'gui' when filename extension is not specified, a format list pops-up
%         options: specific format options, which are usually plot options
%           default is 'view2 axis tight'.
%
% output: f: filename(s) used to save data (char)
% ex:     b=save(a, 'file', 'm');
%         b=save(a, 'file', 'gif', 'axis tight');
%
% Version: $Date$
% See also iFunc, saveas

[filename,format] = saveas(a, varargin{:});

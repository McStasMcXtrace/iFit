function formats = iLoad_ini
% formats = iLoad_ini
%
% definition of import formats to be used by iLoad
%
% each format is specified as a cell with the following elements
% {method, patterns, name, options}
% where:
%   method:   function name to use
%   patterns: list of strings to search in data file. If all found, then method
%             is qualified
%   name:     name of the method/format
%   options:  additional options to pass to the method.
%             If given as a string they are catenated with file name
%             If given as a cell, they are given to the method as additional arguments
%
% formats should be sorted from the most specific to the most general.
% Formats will be tried one after the other, in the given order.
%
% Part of: iFiles utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. June, 2007.

formats= { ...
{ 'looktxt', {'RRRR','AAAA','FFFF','SSSS','IIII'}, 'ILL Data (normal integers)','--headers --fortran --catenate --fast --binary --makerows=IIII --makerows=FFFF'}, ...
{ 'looktxt', {'RRRR','AAAA','FFFF','SSSS','JJJJ'}, 'ILL Data (large integers)', '--headers --fortran --catenate --fast --binary --makerows=JJJJ --makerows=FFFF'}, ...
{ 'looktxt', {'RRRR','AAAA','FFFF','SSSS'},        'ILL (floats only)',         '--headers --fortran --catenate --fast --binary --makerows=FFFF'}, ...
{ 'looktxt', {'SSSS'},                             'ILL Data (general)',        '--headers --fortran --catenate --fast --binary --makerows=FFFF --makerows=JJJJ --makerows=IIII'}, ...
{ 'looktxt', {'ILL TAS data','RRRR','AAAA','VVVV','POLAN'}, 'ILL TAS Data (polarized)' , '--headers --section=PARAM --section=VARIA --section=ZEROS --section=POLAN --metadata=DATA' }, ...
{ 'looktxt', {'ILL TAS data','RRRR','AAAA','VVVV'},'ILL TAS Data' ,  '--headers --section=PARAM --section=VARIA --section=ZEROS --metadata=DATA' }, ...
{ 'looktxt', {}, 'Data (text format with fastest import method)',    '--headers --binary --fast'}, ...
{ 'looktxt', {}, 'Data (text format with fast import method)',       '--headers --binary'}, ...
{ 'looktxt', {}, 'Data (text format)',                               '--headers'}, ...
{ 'importdata',{},'Matlab importer',''}, ...
{ 'wk1read', {}, 'Lotus1-2-3 (first spreadsheet)',''}, ...
{ 'auread',  {}, 'NeXT/SUN (.au) sound',''}, ...
{ 'wavread', {}, 'Microsoft WAVE (.wav) sound',''}, ...
{ 'aviread', {}, 'Audio/Video Interleaved (AVI) ',''}, ...
{ 'cdfread', {}, 'NetCDF',''}, ...
{ 'fitsread',{}, 'FITS',''}, ...
{ 'xlsread', {}, 'Microsoft Excel (first spreadsheet)',''}, ...
{ 'imread',  {}, 'Image',''}, ...
{ 'hdfread', {}, 'HDF4',''}, ...
{ 'hdf5read',{}, 'HDF5',''}, ...
{ 'load',    {}, 'Matlab workspace',''}, ...
{ 'csvread', {}, 'Comma Separated Values',''}, ...
{ 'dlmread', {}, 'Numerical single block',''}, ...
{ 'xmlread', {}, 'XML',''}, ...
{}, ...
};

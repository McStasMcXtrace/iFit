function [filename,format] = saveas(a, varargin)
% f = saveas(s, filename, format) : save iData object into various data formats
%
%   @iData/saveas function to save data sets
%   This function save the content of iData objects. 
%
% input:  s: object or array (iData)
%         filename: name of file to save to. Extension, if missing, is appended (char)
%                   If the filename already exists, the file is overwritten.
%                   If given as filename='gui', a file selector pops-up
%         format: data format to use (char), or determined from file name extension
%           'm'   save as a flat Matlab .m file (a function which returns an iData object or structure)
%           'mat' save as a '.mat' binary file (same as 'save')
%           'hdf' save as an HDF5 data set
%           'nc'  save as NetCDF 
%           'gif','bmp' save as an image (no axes, only for 2D data sets)
%           'png','tiff','jpeg','ps','pdf','ill','eps' save as an image (with axes)
%           'xls' save as an Excel sheet (requires Excel to be installed)
%           'svg' save as Scalable Vector Graphics (SVG) format
%           'vrml' save as Virtual Reality VRML 2.0 file
%           If given as format='gui' and filename extension is not specified, a format list pops-up
%         options: specific format options, which are usually plot options
%           default is 'view2 axis tight'
%
% output: f: filename used to save file(s) (char)
% ex:     b=saveas(a, 'file', 'm');
%         b=saveas(a, 'file', 'svg', 'axis tight');
%
% Contributed code (Matlab Central): 
%   plot2svg:   Juerg Schwizer, 22-Jan-2006 
%
% Version: $Revision: 1.11 $
% See also iData, iData/load, iData/getframe, save

% handle array of objects to save iteratively
if length(a) > 1
  if length(varargin) >= 1, filename_base = varargin{1}; 
  else filename_base = ''; end
  if strcmp(filename_base, 'gui'), filename_base=''; end
  filename = cell(size(a));
  for index=1:length(a(:))
    [filename{index}, format] = saveas(a(index), varargin{:});
    if isempty(filename_base), filename_base = filename{index}; end
    if length(a(:)) > 1
      [path, name, ext] = fileparts(filename_base);
      varargin{1} = [ path name '_' num2str(index) ext ];
    end
  end
  return
end

% default options checks
if nargin < 2, filename = ''; else filename = varargin{1}; end
if isempty(filename), filename = a.Tag; end
if nargin < 3, format=''; else format = varargin{2}; end
if nargin < 4, options=''; else options=varargin{3}; end
if isempty(options) && ndims(a) >= 2, options='view2 axis tight'; end

% supported format list
filterspec = {'*.m',   'Matlab script/function (*.m)'; ...
      '*.mat', 'Matlab binary file (*.mat)'; ...
      '*.pdf', 'Portable Document Format (*.pdf)'; ...
      '*.eps', 'Encapsulated PostScrip (color, *.eps)'; ...
      '*.ps', 'PostScrip (color, *.ps)'; ...
      '*.nc', 'NetCDF (*.nc)'; ...
      '*.hdf', 'Hierarchical Data Format (compressed, *.hdf, *.nx)'; ...
      '*.xls', 'Excel format (requires Excel to be installed, *.xls)'; ...
      '*.csv', 'Comma Separated Values (suitable for Excel, *.csv)'; ...
      '*.png', 'Portable Network Graphics image (*.png)'; ...
      '*.jpg', 'JPEG image (*.jpg)'; ...
      '*.tiff;*.tif', 'TIFF image (*.tif)'; ...
      '*.svg', 'Scalable Vector Graphics (*.svg)'; ...
      '*.wrl', 'Virtual Reality file (*.wrl)'};

% filenape='gui' pops-up a file selector
if strcmp(filename, 'gui')  
  [filename, pathname, filterindex] = uiputfile( ...
       filterspec, ...
        ['Save ' a.Title ' as...'], a.Tag);
  if ~isempty(filename) & filename ~= 0
    ext = filterspec{filterindex,1};
    % check if extension was given
    [f,p,e] = fileparts(filename);
    if isempty(e), filename=[ filename ext(2:end) ]; end
    format=ext(2:end);
  else
    filename=[]; return
  end
end

% format='gui' pops-up a list of available file formats, if not given from file extension
if strcmp(format, 'gui')
  liststring= filterspec{:,2};
  format_index=listdlg('ListString',liststring,'Name',[ 'Select format to save ' filename ], ...
    'PromptString', {'Select format ',['to save file ' filename ]}, ...
    'ListSize', [300 200]);
  if isempty(format_index), return; end
  format = liststring{format_index};
  format = format(3:end);
end

% handle extensions
[path, name, ext] = fileparts(filename);
if isempty(ext) & ~isempty(format), 
  ext = [ '.' format ]; 
  filename = [ filename ext ];
elseif isempty(format) & ~isempty(ext)
  format = ext(2:end);
elseif isempty(format) & isempty(ext) 
  format='m'; filename = [ filename '.m' ];
end

% handle some format aliases
switch format
case 'jpg'
  format='jpeg';
case 'eps'
  format='epsc';
case 'ps'
  format='psc';
case 'netcdf'
  format='nc';
end

% ==============================================================================
% handle specific format actions
switch format
case 'm'  % single m-file Matlab output (text), with the full object description
  str = [ 'function this=' name sprintf('\n') ...
          '% Original data: ' class(a) ' ' inputname(1) ' ' a.Tag sprintf('\n') ...
          '%' sprintf('\n') ...
          '% Matlab ' version ' m-file ' filename ' saved on ' datestr(now) ' with iData/saveas' sprintf('\n') ...
          '% To use import data, type ''' name ''' at the matlab prompt.' sprintf('\n') ...
          '% You will obtain an iData object (if you have iData installed) or a structure.' sprintf('\n') ...
          '%' sprintf('\n') ...
          class2str('this', a) ];
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag ]);
    return
  end
  fprintf(fid, '%s', str);
  fclose(fid);
case 'mat'  % single mat-file Matlab output (binary), with the full object description
  save(filename, a);
case {'hdf','hdf5', 'nc'} % HDF4, HDF5, NetCDF formats: converts fields to double and chars
  [fields, types, dims] = findfield(a);
  towrite={};
  for index=1:length(fields(:)) % get all field names
    val=get(a, fields{index});
    if iscellstr(val), 
      val=val(:);
      val(:, 2)={ ';' }; 
      val=val'; 
      val=[ val{1:(end-1)} ];
    end
    if ~isnumeric(val) & ~ischar(val), continue; end
    if strcmp(format,'nc')
      fields{index} = strrep(fields{index}, '.', '_');
    else
      fields{index} = strrep(fields{index}, '.', filesep);
      if isempty(towrite)
        % initial write wipes out the file
        delete(filename);
        hdf5write(filename, [ filesep 'iData' filesep fields{index} ], val, 'WriteMode', 'overwrite');
      else
        % consecutive calls are appended
        try
          hdf5write(filename, [ filesep 'iData' filesep fields{index} ], val, 'WriteMode', 'append');
        catch
          % object already exists: we skip consecutive write
        end
      end
    end
    if isempty(towrite)
      towrite={ fields{index}, val };
    else
      towrite={ towrite{1:end}, fields{index}, val };
    end
  end
  if strcmp(format,'nc')
    cdfwrite(filename, towrite);
  end
case 'xls'  % Excel file format
  xlswrite(filename, double(a), a.Title);
case 'csv'  % Spreadsheet comma separated values file format
  csvwrite(filename, double(a));
case {'gif','bmp','pbm','pcx','pgm','pnm','ppm','ras','xwd'}  % bitmap images
  if ndims(a) == 2
    a=double(a);
    b=(a-min(a(:)))/(max(a(:))-min(a(:)))*64;
    imwrite(b, jet(64), filename, format);
  else
    f=getframe(a);
    imwrite(f.cdata, jet(64), filename, format);
  end
case 'epsc' % color encapsulated postscript file format, with TIFF preview
  f=figure('visible','off');
  plot(a,options);
  print(f, [ '-depsc -tiff' ], filename);
  close(f);
case {'png','tiff','jpeg','psc','pdf','ill'}  % other bitmap and vector graphics formats (PDF, ...)
  f=figure('visible','off');
  plot(a,options);
  print(f, [ '-d' format ], filename);
  close(f);
case 'fig'  % Matlab figure format
  f=figure('visible','off');
  plot(a,options);
  saveas(f, filename, 'fig');
  close(f);
case 'svg'  % scalable vector graphics format (private function)
  f=figure('visible','off');
  plot(a,options);
  plot2svg(filename, f);
  close(f);
case 'vrml' % VRML format
  f=figure('visible','off');
  plot(a,options);
  vrml(f,filename);
  close(f);
end

% end of iData/saveas

function [filename,format] = saveas(a, varargin)
% f = saveas(s, filename, format, options) : save iData object into various data formats
%
%   @iData/saveas function to save data sets
%     This function saves the content of iData objects. The default format is 'm'.
%   saveas(iData,'formats')
%     prints a list of supported export formats.
%   saveas(iData,'file.ext')            determine file format from the extension
%   saveas(iData,'file','format')       sets file format explicitly
%   saveas(iData,'file','format clean') sets file format explicitly and remove NaN and Inf.
%     To load back a model from an m-file, type its file name at the prompt.
%     To load back a model from an mat-file, type 'load filename.mat' at the prompt.
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
%           'vtk'  save as VTK ascii (<1e5 elements) or binary
%           'x3d'  save as X3D (geometry) file, ascii
%           'xhtml' save as embedded HTML/X3D file (using Flash plugin for rendering)
%           'yaml' save as YAML/JSON format, ascii
%           'xml'  save as an XML file, ascii
%
%           'gui' when filename extension is not specified, a format list pops-up
%         options: specific format options, which are usually plot options
%           default is 'view2 axis tight'
%
% output: f: filename(s) used to save data (char)
% ex:     b=saveas(a, 'file', 'm');
%         b=saveas(a, 'file', 'svg', 'axis tight');
%
% Contributed code (Matlab Central): 
%   plot2svg:   Juerg Schwizer, 22-Jan-2006 
%   iData_private_saveas_hdfnc
%   medf_write
%   fitswrite:  R. G. Abraham, Institute of Astronomy, Cambridge University (1999)
%   stlwrite
%   struct2xml
%   yaml (in Objects)
%
% Version: $Revision$
% See also iData, iData/load, iData/getframe, save

% default options checks
if nargin < 2, filename = ''; else filename = varargin{1}; end
if isempty(filename), filename = [ 'iFit_' a.Tag ]; end
if nargin < 3, format=''; else format = varargin{2}; end
% if the filename is given only as an extension, use it as the format
if nargin == 2 && filename(1) == '.'
  format=filename(2:end);
  filename='';
end

if nargin < 4, options=''; else options=varargin{3}; end
if isempty(options) && any(ndims(a) >= 2), options='view2 axis tight'; end

% supported format list
filterspec = {'*.m',   'Matlab script/function (*.m)'; ...
      '*.dat', 'Flat text file with comments (*.dat)'; ...
      '*.mat', 'Matlab binary file (*.mat)'; ...
      '*.fig', 'Matlab figure (*.fig)'; ...
      '*.pdf', 'Portable Document Format (*.pdf)'; ...
      '*.eps', 'Encapsulated PostScrip (color, *.eps)'; ...
      '*.ps',  'PostScrip (color, *.ps)'; ...
      '*.nc;*.cdf',  'NetCDF (*.nc, *.cdf)'; ...
      '*.hdf;*.hdf5;*.h5','Hierarchical Data Format 5 (*.hdf5, *.h5, *.hdf)'; ...
      '*.hdf4;*.h4', 'Hierarchical Data Format 4 image (*.hdf4)'; ...
      '*.xls', 'Excel format (requires Excel to be installed, *.xls)'; ...
      '*.csv', 'Comma Separated Values (suitable for Excel, *.csv)'; ...
      '*.png', 'Portable Network Graphics image (*.png)'; ...
      '*.jpg', 'JPEG image (*.jpg)'; ...
      '*.tiff;*.tif', 'TIFF image (*.tif)'; ...
      '*.svg', 'Scalable Vector Graphics (*.svg)'; ...
      '*.wrl;*.vrml', 'Virtual Reality file (*.wrl, *.vrml)'; ...
      '*.vtk', 'VTK volume (*.vtk)'; ...
      '*.hdr', 'Analyze volume (*.hdr+img)'; ...
      '*.stl;*.stla;*.stlb', 'Stereolithography geometry (*.stl)'; ...
      '*.off', 'Object File Format geometry (*.off)'; ...
      '*.ply', 'PLY geometry (*.ply)'; ...
      '*.x3d',   'X3D (geometry) file, ascii (*.x3d)'; ...
      '*.xhtml', 'embedded HTML/X3D file (*.html using Flash plugin for rendering)';
      '*.edf', 'EDF ESRF format for 1D and 2D data sets (*.edf)' ; 
      '*.fits;*.fit;*.fts','IAU FITS binary image (*.fits, only for 2D objects)';
      '*.yaml;*.yml;*.json','YAML/JSON interchange format (*.yaml,*.json)';
      '*.xml','XML file (*.xml)' };
if strcmp(filename, 'formats')
  fprintf(1, '       EXT  DESCRIPTION [%s(iData)]\n', mfilename);
  fprintf(1, '-----------------------------------------------------------------\n'); 
  for index=1:size(filterspec,1)
    ext = upper(filterspec{index,1});
    ext = strrep(ext,'.','');
    ext = strrep(ext,'*','');
    fprintf(1,'%10s  %s \n', ext, filterspec{index,2});
  end
  return
end

% filename='gui' pops-up a file selector
if strcmp(filename, 'gui')
  if numel(a) > 1, t=[ num2str(numel(a)) ' objects' ]; else t=get(a,'Title'); end
  [filename, pathname, filterindex] = uiputfile( ...
       filterspec, ...
        ['Save ' t ' as...'], a.Tag);
  if ~isempty(filename) & filename ~= 0
    ext = filterspec{filterindex,1};
    if iscellstr(ext), ext=ext{1}; end
    % check if extension was given
    [f,p,e] = fileparts(filename);
    if isempty(e), 
      filename=[ filename ext(2:end) ];
      format=ext(3:end);
    elseif isempty(format)
      format=e(2:end);
    end
  else
    filename=[]; return
  end
end

% search for option to clean the data set from NaN's and Inf's
index=regexp(format, '\<clean\>');  % search the word
if ~isempty(index)
  index=index(1);
  format(index:(index+length('clean')-1)) = '';
  a = iData_private_cleannaninf(a);
end

% format='gui' pops-up a list of available file formats, if not given from file extension
if any(regexp(format, '\<gui\>'))
  liststring= filterspec(:,2);
  format_index=listdlg('ListString',liststring,'Name',[ 'Select format to save ' filename ], ...
    'PromptString', {'Select format ',['to save file ' filename ]}, ...
    'ListSize', [300 200]);
  if isempty(format_index), return; end
  format = filterspec{format_index,1};
  format = format(3:end);
end

format=lower(strtrim(format));

% handle aliases
switch format
case {'netcdf','nc'}
  format='cdf';
case 'vrml'
  format='wrl';
end

% handle extensions
[path, name, ext] = fileparts(filename);
if isempty(ext) & ~isempty(format), 
  ext = [ '.' format ]; 
  filename = [ filename ext ];
elseif isempty(format) & ~isempty(ext)
  format = ext(2:end);
elseif isempty(format) & isempty(ext) 
  format='mat'; filename = [ filename '.mat' ];
end

% handle array of objects to save iteratively
if numel(a) > 1 & ~strcmp(lower(format),'mat')
  if length(varargin) >= 1, filename_base = varargin{1}; 
  else filename_base = ''; end
  if strcmp(filename_base, 'gui'), filename_base=''; end
  if isempty(filename_base), filename_base='iData'; end
  filename = cell(size(a));
  for index=1:numel(a)
    if numel(a) > 1
      [path, name, ext] = fileparts(filename_base);
      varargin{1} = [ path name '_' num2str(index,'%04d') ext ];
    end
    [filename{index}, format] = saveas(a(index), varargin{:});
  end
  return
end

% handle some format aliases (after extension extraction from file name)
switch format
case 'jpg'
  format='jpeg';
case 'eps'
  format='epsc';
case 'ps'
  format='psc';
case 'netcdf'
  format='cdf';
case 'hdf'
  format='hdf5';
end

% remove NaN values, which are usually not well supported by text based formats

% ==============================================================================
% handle specific format actions
switch format
case 'm'  % single m-file Matlab output (text), with the full object description
  NL = sprintf('\n');
  if ~isdeployed
    str = [ 'function this=' name NL ];
  else
    str = '';
  end
  str = [ str ...
          '% Original data: ' NL ...
          '%   class:    ' class(a) NL ...
          '%   variable: ' inputname(1) NL ...
          '%   tag:      ' a.Tag NL ...
          '%   label:    ' a.Label NL ...
          '%   source:   ' a.Source NL ... 
          '%' NL ...
          '% Matlab ' version ' m-file ' filename ' saved on ' datestr(now) ' with iData/saveas' NL ...
          '% To use/import data, type ''' name ''' at the matlab prompt.' NL ...
          '% You will obtain an iData object (if you have iData installed) or a structure.' NL ...
          '%' NL ...
          class2str('this', a) ];
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag 'in format ' format ]);
    disp(message)
    return
  end
  fprintf(fid, '%s', str);
  fclose(fid);
  if isdeployed
    disp([ 'Warning: The standalone/deployed version of iFit does not allow to read back' NL ...
           '  function definitions. This m-file has been converted to a script that you can' NL ...
           '  import as "this" by typing: run ' filename ]);
  end
case 'dat'  % flat text file with commented blocks, in the style of McStas/PGPLOT
  NL = sprintf('\n');
  str = [ '# Format: data with text headers' NL ...
          '# URL: ifit.mccode.org' NL ...
          '# Creator: iFit/@iData/saveas - ' version(a) NL ...
          '# Title: ' a.Title NL ...
          '# Label: ' a.Label NL ...
          '# DisplayName: ' a.DisplayName NL ...
          '# User: ' a.User NL ...
          '# CreationDate: ' get(a,'Date') NL ...
          '# ModificationDate: ' get(a,'ModificationDate') NL ...
          '# Tag: ' a.Tag NL ...
          '# ' NL ...
          class2str('', a, 'flat') ];
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag 'in format ' format ]);
    disp(message)
    return
  end
  fprintf(fid, '%s', str);
  fclose(fid);
case 'mat'  % single mat-file Matlab output (binary), with the full object description
  if ~isempty(inputname(1))
    eval([ inputname(1) '= a;' ]);
    save(filename, inputname(1));
  else
    eval([ a.Tag '= a;' ]);
    save(filename, a.Tag);
  end
case {'hdf5', 'nc','cdf', 'nx','h5'} % HDF4, HDF5, NetCDF formats: converts fields to double and chars
  filename = iData_private_saveas_hdfnc(a, filename, format); % inline function (below)
case 'edf'  % EDF ESRF format
  filename = medfwrite(a, filename); % in private
case 'vtk'  % VTK volume
  filename = iData_private_saveas_vtk(a, filename);
case 'hdr'  % Analyze volume
  filename = iData_private_saveas_analyze(a, filename);
case {'fits','fit','fts'} % FITS image
  if ndims(a) == 2
    a = double(a);
    fitswrite(a, filename);
  end
case 'xls'  % Excel file format
  xlswrite(filename, double(a), a.Title);
case 'csv'  % Spreadsheet comma separated values file format
  csvwrite(filename, double(a));
case {'gif','bmp','pbm','pcx','pgm','pnm','ppm','ras','xwd','hdf4'}  % bitmap images
  if ndims(a) == 2 
    b=getaxis(a,0); % Signal/Monitor
    if abs(log10(size(b,1)) - log10(size(b,2))) > 1
      x = round(linspace(1, size(b,1), max(size(b,1), 1024)));
      y = round(linspace(1, size(b,2), max(size(b,2), 1024)));
      b = b(x,y);
    end
    b=(b-min(b(:)))/(max(b(:))-min(b(:)))*64;
  else
    f=getframe(a);
    b = f.cdata;
  end
  switch format
  case 'gif'
    imwrite(b, jet(64), filename, format, 'Comment',char(a));
  otherwise
    imwrite(b, jet(64), filename, format);
  end
case 'epsc' % color encapsulated postscript file format, with TIFF preview
  f=figure('visible','off');
  plot(a,options);
  print(f, '-depsc', '-tiff', filename);
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
case {'vrml','wrl'} % VRML format
  f=figure('visible','off');
  h = plot(a,options);
  g = gca;
  vrml(g,filename);
  close(f);
case {'x3d','xhtml'} % X3D/XHTML format
  f=figure('visible','off');
  h = plot(a,options);
  titl = strrep(char(a),'<','[');
  titl = strrep(titl,'>',']');
  desc = evalc('disp(a)');
  desc = strrep(desc,'<','[');
  desc = strrep(desc,'>',']');
  figure2xhtml(filename, f, struct('interactive',true,'output', format,'title',titl,'Description',desc));
  close(f);
case {'stl','stla','stlb','off','ply'} % STL ascii, binary, PLY, OFF
  if ndims(a) == 1    iData_private_warning(mfilename,[ 'Object ' inputname(1) ' ' a.Tag ' does not seem to be exportatble as a ' format ' file. Ignoring.' ]);
    return
  else
    if ndims(a) == 2
      [x, xlab] = getaxis(a,2); x=double(x);
      [y, ylab] = getaxis(a,1); y=double(y);
      [z, clab] = getaxis(a,0); c=double(c);
    elseif ndims(a) >= 3
      [x, xlab] = getaxis(a,2); x=double(x);
      [y, ylab] = getaxis(a,1); y=double(y);
      [z, zlab] = getaxis(a,3); z=double(z);
    end
    if any(strcmp(format, {'stl','stlb'}))
      mode = 'binary';
    else
      mode = 'ascii';
    end
    
    % get Title
    T   = a.Title; if ~ischar(T), T=char(T); end
    if ~isvector(T), T=transpose(T); T=T(:)'; end
    T   = regexprep(T,'\s+',' '); % remove duplicated spaces
    if length(T) > 69, T=[ T(1:60) '...' T((end-8):end) ]; end
    % get the faces and vertices
    v = [x(:) y(:) z(:)];
    f = delaunay(x(:),y(:));
    %f=delaunayn(v,{'Qx','Qv','Tv', 'Qt','Qbb','Qc','Qz'});
    [v, indexm, indexn] =  unique(v, 'rows');
    f = indexn(f);
    if strncmp(format,'stl',3)  % STL format
      stlwrite(filename, f,v, 'Mode', mode, 'Title', T);
    else                        % OFF and PLY formats
      [fid, message]=fopen(filename,'w+');
      if fid == -1
        iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag 'in format ' format ]);
        disp(message)
        return
      end
      % write header
      if strcmp(format,'ply')   % OFF format
        fprintf(fid,'ply\nformat ascii 1.0\ncomment This is a PLY file format for %s\nelement vertex %d\n', ...
          T, size(v,1));
        fprintf(fid,'property float x\nproperty float y\nproperty float z\nelement face %d\n', size(f,1));
        fprintf(fid,'property list uchar int vertex_indices\nend_header\n');
      else                      % PLY
        fprintf(fid,'OFF\n%i %i 0\n',...
          size(v,1), size(f,1));
      end
      % write data
      str = num2str(v,5); str(:,end+1) = sprintf('\n');
      str = str'; str = str(:)';
      fprintf(fid,'%s', str);
      f = [ (size(f,2)+1)*ones(size(f,1),1) f ];
      str = num2str(f-1); str(:,end+1) = sprintf('\n');
      str = str'; str = str(:)';
      fprintf(fid,'%s', str);
      if strcmp(format,'off')
        fprintf(fid,'# This is an Object File Format (geomview) for %s\n', T);
       end
      fclose(fid);
    end
  end
case {'yaml','yml','json'}
  YAML.write( filename, struct(a) ); % YAML object is in iFit/Objects
case {'xml'}
  struct2xml(struct(a), filename);   % in private
otherwise
  iData_private_warning(mfilename,[ 'Export of object ' inputname(1) ' ' a.Tag ' into format ' format ' is not supported. Ignoring.' ]);
  filename = [];
end

% end of iData/saveas



function [filename,format] = saveas(a, filename, format, options)
% f = saveas(s, filename, format, options) : save iData object into various data formats
%
%   @iData/saveas function to save data sets
%     This function saves the content of iData objects. The default format is 'm'.
%   saveas(iData,'formats')
%     prints a list of supported export formats.
%   saveas(s,'file.ext')            determine file format from the extension
%   saveas(s,'file','format')       set file format explicitly
%   saveas(s,'file','format clean') set file format explicitly and remove NaN and Inf.
%   saveas(s,'file','format data')  save only the 'Data' part of the object. 
%   saveas(s, get(s,'Source'), 'format')  save with the same base filename(s)
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
%           tight, hide_axes, interp, view2, view3, transparent, light, clabel, colorbar, whole
%           default is 'view2 axis tight'
%           For XHTML export, the additional argument can be a string with the isosurface
%           level and options 'axes' to display axes, 'auto' to rescale as a cube.
%
% output: f: filename(s) used to save data (char)
% ex:     b=saveas(a, 'file', 'm');
%         b=saveas(a, 'file', 'svg', 'axis tight');
%         b=saveas(a, 'file', 'hdf data');
%
% Version: $Date$
% See also iData, iData/load, iData/getframe, save, iData/plot

% Contributed code (Matlab Central): 
%   plot2svg:   Juerg Schwizer, 22-Jan-2006 
%   medf_write
%   fitswrite:  R. G. Abraham, Institute of Astronomy, Cambridge University (1999)
%   stlwrite
%   struct2xml
%   yaml (in Objects)
%   mat2json
%   exportToPPTX
%   writeNPY
%   mesh2kml
%
%   iData_private_saveas_hdfnc

% default options checks
if nargin < 2, filename = ''; end
if nargin < 3, format='';     end
% if the filename is given only as an extension, use it as the format
if nargin == 2 && filename(1) == '.'
  format=filename(2:end);
  filename='';
end

if nargin < 4, options=''; end
if isempty(options) && any(ndims(a) >= 2), options='view2 axis tight'; end

% supported format list
filterspec = { ...
      '*.avi', 'Audio Video Interleave (AVI) movie (*.avi)'; ...
      '*.cdf',  'CDF (*.cdf)'; ...
      '*.csv', 'Comma Separated Values (suitable for Excel, *.csv)'; ...
      '*.dae', 'Collada model (*.dae)'; ...
      '*.dat', 'Flat text file with comments (*.dat)'; ...
      '*.edf', 'EDF ESRF format for 1D and 2D data sets (*.edf)' ; 
      '*.eps', 'Encapsulated PostScript (color, *.eps)'; ...
      '*.fig', 'Matlab figure (*.fig)'; ...
      '*.fits;*.fit;*.fts','IAU FITS binary image (*.fits, only for 2D objects)';
      '*.hdf;*.hdf5;*.h5;*.nxs;*.n5','Hierarchical Data Format 5 (*.hdf5, *.h5, *.hdf)'; ...
      '*.hdf4;*.h4;*.nxs;*.n4', 'Hierarchical Data Format 4 image (*.hdf4)'; ...
      '*.hdr', 'Analyze volume (*.hdr+img)'; ...
      '*.html;*.htm','Hypertext Markup Language document (*.html)'; ...
      '*.inx;*.ino','ILL Inelastic Neutron Scattering data (*.inx, *.ino)'; ...
      '*.jpg;*.jpeg', 'JPEG image (*.jpg)'; ...
      '*.json', 'JSON JavaScript Object Notation (*.json)'; ...
      '*.kml', 'GoogleEarth model (*.kml)'; ...
      '*.m',   'Matlab script/function (*.m)'; ...
      '*.mat', 'Matlab binary file (*.mat, serialized)'; ...
      '*.mrc', 'MRC map file (*.mrc)'; ...
      '*.nc;*.cdf',  'NetCDF (*.nc, *.cdf)'; ...
      '*.nii', 'NiFti volume (*.nii)'; ...
      '*.npy', 'Numpy array NPY (*.npy)'; ...
      '*.off', 'Object File Format geometry (*.off)'; ...
      '*.ply', 'PLY geometry (*.ply)'; ...
      '*.pdf', 'Portable Document Format (*.pdf)'; ...
      '*.ps',  'PostScript (color, *.ps)'; ...
      '*.png', 'Portable Network Graphics image (*.png)'; ...
      '*.ppm;*.pbm;*.pgm;*.pnm','Portable Anymap format (*.pnm)'; ...
      '*.spe', 'ISIS SPE Mslice/Horace (*.spe)'; ...
      '*.sqw', 'McStas SQW Isotropic S(q,w) (*.sqw)'; ...
      '*.stl;*.stla;*.stlb', 'Stereolithography geometry (*.stl)'; ...
      '*.svg', 'Scalable Vector Graphics (*.svg)'; ...
      '*.tiff;*.tif', 'TIFF image (*.tif)'; ...
      '*.vtk', 'VTK volume (*.vtk)'; ...
      '*.wrl;*.vrml', 'Virtual Reality file (*.wrl, *.vrml)'; ...
      '*.x3d',   'X3D (geometry) file, ascii (*.x3d)'; ...
      '*.xhtml', 'embedded HTML/X3D file (*.html using Flash plugin for rendering)'; ...
      '*.xls', 'Excel format (requires Excel to be installed, *.xls)'; ...
      '*.xml','XML file (*.xml)'; ...
      '*.yaml;*.yml','YAML interchange format (*.yaml)' };
if strcmp(filename, 'formats')
  fprintf(1, '       EXT  DESCRIPTION [%s(iData)]\n', mfilename);
  fprintf(1, '-----------------------------------------------------------------\n'); 
  for index=1:size(filterspec,1)
    ext = upper(filterspec{index,1});
    ext = strrep(ext,'.','');
    ext = strrep(ext,'*','');
    fprintf(1,'%10s  %s \n', ext, filterspec{index,2});
  end
  filename = filterspec;
  return
end

% filename='gui' pops-up a file selector
if strcmp(filename, 'gui')
  if numel(a) > 1, t=[ num2str(numel(a)) ' objects' ]; else t=get(a,'Title'); end
  [filename, pathname, filterindex] = uiputfile( ...
       filterspec, ...
        ['Save ' t ' as...'], a.Tag);
  if ~isempty(filename) & filename ~= 0
    ext = strtok(filterspec{filterindex,1},' *;');
    if iscell(ext) && ischar(ext{1}), ext=ext{1}; end
    % check if extension was given
    [f,p,e] = fileparts(filename);
    if isempty(e), 
      filename=[ filename ext ];
      format=ext;
    elseif isempty(format)
      format=e;
    end
    filename = strcat(pathname, filesep, filename);
  else
    filename=[]; return
  end
end

if isempty(format) && ~isempty(filename) && any(~cellfun(@isempty,strfind(filterspec(:,1),strtok(filename, ' ;*'))))
  format = filename; filename = '';
end
if ~isempty(format) && format(1) == '.', format=format(2:end); end
if strcmp(format, 'mantid')
  format = 'nxs mantid data';
end
if strcmp(format, 'lamp')
  format = 'nxs lamp data';
end

% search for option to clean the data set from NaN's and Inf's
index=regexp(format, '\<clean\>');  % search the word
if ~isempty(index)
  a = iData_private_cleannaninf(a);
end

% convert data set as Mantid Processed Workspace when requested
index=regexp(format, '\<mantid\>');  % search the word
if ~isempty(index)
  a = iData_private_2mantid(a);
end

% convert data set as LAMP Processed Workspace when requested
index=regexp(format, '\<lamp\>');  % search the word
if ~isempty(index)
  a = iData_private_2lamp(a);
end

% search the word 'data' to only save object Data property (for HDF,CDF,NetCDF)
index=regexp(format, '\<data\>');  % search the word
if ~isempty(index)
  root   = 'Data';
else root='';
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
case 'netcdf'
  format='nc';
case 'vrml'
  format='wrl';
case 'mantid'
  format='hdf5 mantid data';
end

% handle extensions
if ~ischar(filename), filename=char(filename); end

if ~isdir(filename)
  [Path, name, ext] = fileparts(filename);
  file_isdir = '';
else  % isdir
  name = ''; ext=''; 
  file_isdir = filename; filename = '';
end
if isempty(ext) && ~isempty(format), 
  ext = [ '.' strtok(format, ' ;*') ]; 
  filename = [ filename ext ];
elseif isempty(format) && ~isempty(ext)
  format = ext(2:end);
elseif isempty(format) && isempty(ext) 
  format='mat'; filename = [ filename '.mat' ];
end

if isempty(filename) || isempty(name), 
  filename = [ 'iFit_' a(1).Tag '.' strtok(format, ' ;*') ]; 
  if ~isempty(file_isdir), filename = fullfile(file_isdir, filename); end
  name = filename; 
end


formatShort = strtok(format, ' ;*.');
% handle array of objects to save iteratively, except for file formats that support
% multiple entries: HTML MAT
if numel(a) > 1
  filename_base = filename;
  if iscellstr(filename) && numel(filename) == numel(a)
    for index=1:numel(a)
      [Path, name, ext] = fileparts(filename{index});
      [filename{index}, format] = saveas(a(index), [name '.' strtok(format)], format, options);
    end
  else
    if strcmp(filename_base, 'gui'), filename_base=''; end
    if isempty(filename_base),       filename_base='iFit_'; end
    filename = cell(size(a));
    for index=1:numel(a)
      if ~strcmpi(formatShort, 'html') && ~strcmpi(formatShort, 'mat')
        [Path, name, ext] = fileparts(filename_base);
        this_filename = fullfile(Path, [ name '_' num2str(index,'%04d') ext ]);
      else
        if index == 1 && ~isempty(dir(filename_base))
          delete(filename_base);
        end
        this_filename = filename_base;
        format = [ format ' ' root ];
      end
      [filename{index}, format] = saveas(a(index), this_filename, format, options);
    end
  end
  return
end


% handle some format aliases (after extension extraction from file name)
switch formatShort
case 'jpg'
  formatShort='jpeg';
case 'eps'
  formatShort='epsc';
case 'ps'
  formatShort='psc';
case 'netcdf'
  formatShort='nc';
end

% remove NaN values, which are usually not well supported by text based formats

% ==============================================================================
% handle specific format actions
try
  switch formatShort
  case 'avi'
    filename = iData_private_saveas_avi(a, filename, options);
  case 'csv'  % Spreadsheet comma separated values file format
    csvwrite(filename, double(a));
  case 'dae'  % Collada DAE
    if ndims(a) == 2
      x = getaxis(a,1);
      y = getaxis(a,2);
      z = getaxis(a,0);
      mesh2kml(x,y,z,filename, ...
        'color','red','alpha',0.7);
    end
  case 'dat'  % flat text file with commented blocks, in the style of McStas/PGPLOT
    filename = iData_private_saveas_dat(a, filename);
  case 'edf'  % EDF ESRF format
    filename = medfwrite(a, filename); % in private
  case 'epsc' % color encapsulated postscript file format, with TIFF preview
    f=figure('visible','off');
    plot(a,options);
    set(f,'renderer','zbuffer');
    print(f, '-depsc', filename); % tiff preview may be broken on some GPU
    close(f);
  case 'fig'  % Matlab figure format
    f=figure('visible','off');
    plot(a,options);
    saveas(f, filename, 'fig');
    close(f);
  case {'fits','fit','fts'} % FITS image
    if ndims(a) == 2
      b = real(double(iData_private_cleannaninf(a,1))); % Signal/Monitor
      fitswrite(b', filename);
    else
      disp([ mfilename ': Export into ' format ' is only possible for 2D objects, not for ' num2str(ndims(a)) 'D. Use resize to change dimensionality. Ignoring.' ]) 
    end
  case {'gif','bmp','pbm','pcx','pgm','pnm','ppm','ras','xwd','hdf4','tiff','png','art'}  % bitmap images
    if ndims(a) == 2 && ~isempty(root) % 'data' option
      b = real(double(iData_private_cleannaninf(a,1)));
      b = round((b-min(b(:)))/(max(b(:))-min(b(:)))*256);
      b = flipud(b);
    else
      f = getframe(a,[],options);
      b = f.cdata;
      if  strcmp(formatShort(1:3),'jpe')
          b=sum(b,3);
      end
    end
    if strcmp(formatShort,'hdf4'), formatShort='hdf'; end
    if ~isempty(b)
      switch formatShort
      case {'png'}
        try
          imwrite(b, jet(256), filename, formatShort, 'Comment',char(a),'Mode','lossless', ...
            'Source',a.Source, 'Software', version(iData), 'CreationTime', a.Date);
        catch
          imwrite(b, jet(256), filename, formatShort, 'Comment',char(a),'Mode','lossless', ...
            'Source',a.Source, 'Software', version(iData));
        end
      case {'gif'}
        imwrite(b, filename, formatShort, 'Comment',char(a));
      case 'tiff'
        imwrite(b, jet(256), filename, formatShort, 'Description',char(a));
      case 'art'
        textart(b, filename); % in private
      otherwise
        try
          imwrite(b, jet(256), filename, formatShort);
        catch
          imwrite(b, filename, formatShort);
        end
      end
    else
      % rendering with getframe/imwrite failed. Fall back to saveas.
      f = figure('visible','off');
      b = plot(a, options);
      set(f,'Renderer','zbuffer')
      saveas(f, filename, formatShort);
      close(f);
    end
  case {'hdf','hdf5','h5','nx','nxs','n5','nc','cdf'} % HDF5, CDF, NetCDF formats: converts fields to double and chars
    filename = iData_private_saveas_hdfnc(a, filename, formatShort, root); % private function
  case 'hdr'  % Analyze volume
    filename = iData_private_saveas_analyze(a, filename);
  case {'html','htm'}
    % create a folder with the HTML doc, figures
    filename = iData_private_saveas_html(a, filename, [ format ' ' root ' ' options ]);
  case {'inx','mcstas','sqw','spe'}  % ILL INX / McStas / SPE
      a = iData_Sqw2D(a);
      filename = saveas(a, filename, formatShort);
  case 'json'
    mat2json(struct(a), filename );    % in private
  case 'kml'  % Google KML
    if ndims(a) == 2
      x = getaxis(a,1);
      y = getaxis(a,2);
      z = getaxis(a,0);
      mesh2kml(x,y,z,filename, ...
        'color','red','position',[-33.85622 151.21535 10],'alpha',0.7);
    end
  case 'm'  % single m-file Matlab output (text), with the full object description
    filename = iData_private_saveas_m(a, filename, name, options);
  case 'mat'  % single mat-file Matlab output (binary), with the full object description
    % serialize for much faster save
    a.Data = hlp_serialize(a.Data);
    varg = { filename };
    if ~isempty(inputname(1))
      eval([ inputname(1) '= a;' ]);
      varg{2} = inputname(1);
    else
      eval([ a.Tag '= a;' ]);
      varg{2} = a.Tag;
    end
    if isempty(dir(filename))
      disp([ mfilename ': The file ' filename ' has been serialized. You MUST import it with load(iData, ''' filename ''')' ])
    else varg{3} = '-append';
    end
    save(varg{:});
  case 'mrc'  % MRC map file
    WriteMRC(getaxis(iData_private_cleannaninf(a),0),1,filename);  % in private
  case 'nii'  % NifTi volume
    filename = iData_private_saveas_nii(a, filename);
  case 'npy'  % Numpy binary array NPY
    writeNPY(getaxis(iData_private_cleannaninf(a),0),filename);
  case {'ppt','pptx'}
    f = figure('visible','off');
    set(f,'renderer','zbuffer');
    b = plot(a, options);
    exportToPPTX('new');
    exportToPPTX('addslide');
    exportToPPTX('addtext', char(a));
    exportToPPTX('addslide');
    exportToPPTX('addpicture',f,'Scale','maxfixed');
    close(f);
    exportToPPTX('addnote',char(a));
    exportToPPTX('save',filename);
    exportToPPTX('close');
  case {'psc','pdf','ill','jpeg'}  % other bitmap and vector graphics formats (PDF, ...)
    f=figure('visible','off');
    plot(a,options);
    set(f,'Renderer','zbuffer')
    saveas(f, filename, formatShort);
    close(f);
  case {'stl','stla','stlb','off','ply'} % STL ascii, binary, PLY, OFF
    filename = iData_private_saveas_stl(a, filename, formatShort);
  case 'svg'  % scalable vector graphics format (private function)
    f=figure('visible','off');
    plot(a,options);
    set(f,'Renderer','zbuffer')
    try
      saveas(f, filename, 'svg');
    catch
      plot2svg(filename, f);
    end
    close(f);
  case 'vtk'  % VTK volume
    filename = iData_private_saveas_vtk(a, filename);
  case {'vrml','wrl'} % VRML format
    f=figure('visible','off');
    h = plot(a,options);
    set(f,'renderer','zbuffer');
    g = gca;
    vrml(g,filename);
    close(f);
  case {'x3d','xhtml'} % X3D/XHTML format
    filename = iData_private_saveas_x3d(a, filename, formatShort, options);
  case 'xls'  % Excel file format
    xlswrite(filename, double(a), a.Title);
  case {'xml'}
    struct2xml(struct(a), filename);   % in private
  case {'yaml','yml'}
    if usejava('jvm')
      YAML.write( filename, struct(a) ); % YAML object is in iFit/Objects
    end
  otherwise
    iData_private_warning(mfilename,[ 'Export of object ' inputname(1) ' ' a.Tag ' into format ' format ' is not supported. Ignoring.' ]);
    filename = [];
  end
catch ME
  disp(getReport(ME))
  iData_private_warning(mfilename,[ 'Export of object ' inputname(1) ' ' a.Tag ' into format ' format ' failed. Ignoring.' ]);
  filename = [];
end

% end of iData/saveas



function [filename,format] = saveas(a, filename, format, options)
% f = saveas(s, filename, format, options) : save iData object into various data formats
%
%   @iData/saveas function to save data sets
%     This function saves the content of iData objects. The default format is 'm'.
%   saveas(iData,'formats')
%     prints a list of supported export formats.
%   saveas(iData,'file.ext')            determine file format from the extension
%   saveas(iData,'file','format')       set file format explicitly
%   saveas(iData,'file','format clean') set file format explicitly and remove NaN and Inf.
%   saveas(iData,'file','format data')  save only the 'Data' part of the object. 
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
%           'art'  save as ASCII art
%           'cdf'  save as CDF (not recommended)
%           'hdf5' save as an HDF5 data set ('nxs','n5','h5' also work)
%           'm'    save as a flat Matlab .m file (a function which returns an iData object or structure)
%           'mantid' save as Mantid Processed Workspace, i.e. 'nxs mantid data' (HDF5)
%           'lamp' save as LAMP Processed Workspace, i.e. 'nxs lamp data'
%           'mat'  save as a serialized '.mat' binary file (fast 'save', DEFAULT)
%           'nc'   save as NetCDF
%         as well as other lossy formats
%           'csv'  save as a comma separated value file
%           'dat'  save as Flat text file with comments
%           'edf'  EDF ESRF format for 1D and 2D data sets
%           'fig'  save as a Matlab figure
%           'fits' save as FITS binary image (only for 2D objects)
%           'gif','bmp','png','tiff','jpeg' save as an image (no axes, only for 2D data sets)
%           'hdf4' save as an HDF4 image
%           'hdr'  save as HDR/IMG Analyze MRI volume (3D/4D)
%           'html' save as Hypertext Markup Language document
%           'json' save as JSON JavaScript Object Notation, ascii
%           'mrc'  save as MRC map file (3/4D)
%           'nii'  save as NifTi Neuroimaging Informatics Technology Initiative (3/4D)
%           'ps','pdf','ill','eps' save as an image (with axes)
%           'ppm','pgm','pbm'
%           'off'  save as Object File Format (geometry), ascii
%           'ply'  save as PLY (geometry), ascii
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
if isempty(filename), filename = [ 'iFit_' a(1).Tag ]; end

if nargin < 4, options=''; end
if isempty(options) && any(ndims(a) >= 2), options='view2 axis tight'; end

% supported format list
filterspec = { ...
      '*.csv', 'Comma Separated Values (suitable for Excel, *.csv)'; ...
      '*.dat', 'Flat text file with comments (*.dat)'; ...
      '*.edf', 'EDF ESRF format for 1D and 2D data sets (*.edf)' ; 
      '*.eps', 'Encapsulated PostScrip (color, *.eps)'; ...
      '*.fig', 'Matlab figure (*.fig)'; ...
      '*.fits;*.fit;*.fts','IAU FITS binary image (*.fits, only for 2D objects)';
      '*.hdf;*.hdf5;*.h5;*.nxs;*.n5','Hierarchical Data Format 5 (*.hdf5, *.h5, *.hdf)'; ...
      '*.hdf4;*.h4;*.nxs;*.n4', 'Hierarchical Data Format 4 image (*.hdf4)'; ...
      '*.hdr', 'Analyze volume (*.hdr+img)'; ...
      '*.html;*.htm','Hypertext Markup Language document (*.html)'; ...
      '*.jpg;*.jpeg', 'JPEG image (*.jpg)'; ...
      '*.json', 'JSON JavaScript Object Notation (*.json)'; ...
      '*.m',   'Matlab script/function (*.m)'; ...
      '*.mat', 'Matlab binary file (*.mat, serialized)'; ...
      '*.mrc', 'MRC map file (*.mrc)'; ...
      '*.nc;*.cdf',  'NetCDF (*.nc, *.cdf)'; ...
      '*.nii','NiFti volume (*.nii)'; ...
      '*.cdf',  'CDF (*.cdf)'; ...
      '*.off', 'Object File Format geometry (*.off)'; ...
      '*.ply', 'PLY geometry (*.ply)'; ...
      '*.pdf', 'Portable Document Format (*.pdf)'; ...
      '*.ps',  'PostScrip (color, *.ps)'; ...
      '*.png', 'Portable Network Graphics image (*.png)'; ...
      '*.ppm;*.pbm;*.pgm;*.pnm','Portable Anymap format (*.pnm)'; ...
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
    if iscell(ext) && ischar(ext{1}), ext=ext{1}; end
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

if isempty(format) && ~isempty(filename) && any(~cellfun(@isempty,strfind(filterspec(:,1),strtok(filename))))
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
  index=index(1);
  format(index:(index+length('clean')-1)) = '';
  a = iData_private_cleannaninf(a);
end

% convert data set as Mantid Processed Workspace when requested
index=regexp(format, '\<mantid\>');  % search the word
if ~isempty(index)
  index=index(1);
  format(index:(index+length('mantid')-1)) = '';
  a = iData_private_2mantid(a);
end

% convert data set as LAMP Processed Workspace when requested
index=regexp(format, '\<lamp\>');  % search the word
if ~isempty(index)
  index=index(1);
  format(index:(index+length('lamp')-1)) = '';
  a = iData_private_2lamp(a);
end

% search the word 'data' to only save object Data property (for HDF,CDF,NetCDF)
index=regexp(format, '\<data\>');  % search the word
if ~isempty(index)
  index=index(1);
  format(index:(index+length('data')-1)) = '';
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
[Path, name, ext] = fileparts(filename);
if isempty(ext) && ~isempty(format), 
  ext = [ '.' format ]; 
  filename = [ filename ext ];
elseif isempty(format) && ~isempty(ext)
  format = ext(2:end); filename = '';
elseif isempty(format) && isempty(ext) 
  format='mat'; filename = [ filename '.mat' ];
end

% handle array of objects to save iteratively
if numel(a) > 1 && ~any(strcmp(lower(format),'mat'))
  filename_base = filename;
  if strcmp(filename_base, 'gui'), filename_base=''; end
  if isempty(filename_base),       filename_base='iFit_'; end
  filename = cell(size(a));
  for index=1:numel(a)
    if numel(a) > 1
      [path, name, ext] = fileparts(filename_base);
      this_filename = [ path name '_' num2str(index,'%04d') ext ];
    end
    [filename{index}, format] = saveas(a(index), this_filename, format, options);
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
  format='nc';
end

% remove NaN values, which are usually not well supported by text based formats

% ==============================================================================
% handle specific format actions
switch strtok(format)
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
          class2str('this', a, options) ];
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
            '# URL: http://ifit.mccode.org' NL ...
            '# Creator: iFit/@iData/saveas - ' version(a) NL ...
            '# Title: ' a.Title NL ...
            '# Label: ' a.Label NL ...
            '# DisplayName: ' a.DisplayName NL ...
            '# User: ' a.User NL ...
            '# CreationDate: ' get(a,'Date') NL ...
            '# ModificationDate: ' get(a,'ModificationDate') NL ...
            '# Tag: ' a.Tag NL ];
  n = sprintf('# Type: %dD', ndims(a));
  if ~isvector(a) % histogram
    n = [ n ' data set ' mat2str(size(a)) ];
    str = [ str n NL class2str('', a, 'flat') ];
  else            % event list
    dat = zeros(prod(size(a)), ndims(a)+1);
    lab = '# Data:';
    for index=1:ndims(a)
      dat(:,index) = getaxis(a, index); % store axes first
      this_lab = label(a, index);
      if isempty(this_lab), this_lab=getaxis(a, num2str(index)); end
      lab = [ lab ' ' this_lab ];
    end
    dat(:, ndims(a)+1) = getaxis(a, 0);
    lab = [ lab ' ' 'Signal' ];
    n = [ n ' event list ' mat2str(size(a)) ];
    str = [ str n NL lab NL class2str('', dat, 'flat') ];
  end
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag 'in format ' format ]);
    disp(message)
    return
  end
  fprintf(fid, '%s', str);
  fclose(fid);
case 'mat'  % single mat-file Matlab output (binary), with the full object description
  % serialize for much faster save
  a.Data = hlp_serialize(a.Data);
  if ~isempty(inputname(1))
    eval([ inputname(1) '= a;' ]);
    save(filename, inputname(1));
  else
    eval([ a.Tag '= a;' ]);
    save(filename, a.Tag);
  end
case {'hdf','hdf5','h5','nx','nxs','n5','nc','cdf'} % HDF5, CDF, NetCDF formats: converts fields to double and chars
  filename = iData_private_saveas_hdfnc(a, filename, format, root); % private function
case 'edf'  % EDF ESRF format
  filename = medfwrite(a, filename); % in private
case 'vtk'  % VTK volume
  filename = iData_private_saveas_vtk(a, filename);
case 'nii'  % NifTi volume
  filename = iData_private_saveas_nii(a, filename);
case 'hdr'  % Analyze volume
  filename = iData_private_saveas_analyze(a, filename);
case 'mrc'  % MRC map file
  WriteMRC(getaxis(a,0),1,filename);  % in private
case {'fits','fit','fts'} % FITS image
  if ndims(a) == 2
    a = double(a);
    fitswrite(a, filename);
  else
    disp([ mfilename ': Export into ' format ' is only possible for 2D objects, not for ' num2str(ndims(a)) 'D. Use resize to change dimensionality. Ignoring.' ]) 
  end
case 'xls'  % Excel file format
  xlswrite(filename, double(a), a.Title);
case 'csv'  % Spreadsheet comma separated values file format
  csvwrite(filename, double(a));
case {'gif','bmp','pbm','pcx','pgm','pnm','ppm','ras','xwd','hdf4','tiff','png','art'}  % bitmap images
  if ndims(a) == 2 && (any(regexp(options, '\<hide_axes\>')))
    b=getaxis(a,0); % Signal/Monitor
    b=round((b-min(b(:)))/(max(b(:))-min(b(:)))*256);
  else
    f = getframe(a,[],options);
    b = f.cdata;
    if  strcmp(format(1:3),'jpe')
        b=sum(b,3);
    end
  end
  switch format
  case {'png'}
    imwrite(b, jet(256), filename, format, 'Comment',char(a),'Mode','lossless');
  case 'tiff'
    imwrite(b, jet(256), filename, format, 'Description',char(a));
  case 'art'
    textart(b, filename); % in private
  otherwise
    if strcmp(format,'hdf4'), format='hdf'; end
    imwrite(b, jet(256), filename, format);
  end
case 'epsc' % color encapsulated postscript file format, with TIFF preview
  f=figure('visible','off');
  plot(a,options);
  print(f, '-depsc', '-tiff', filename);
  close(f);
case {'psc','pdf','ill','jpeg'}  % other bitmap and vector graphics formats (PDF, ...)
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
  titl = char(a);
  titl(titl=='<')='[';
  titl(titl=='>')=']';
  desc = evalc('disp(a)');
  desc(desc=='<')='[';
  desc(desc=='>')=']';
  if ndims(a) <= 2
    f=figure('visible','off');
    if ischar(options)
      h = plot(reducevolume(a),options); % make sure file is not too big
    else
      h = plot(reducevolume(a));
    end
    figure2xhtml(filename, f, struct('interactive',true, ...
      'output', format,'title',titl,'Description',desc));
    close(f);
  else
    [x, xlab] = getaxis(a,2); x=double(x);
    [y, ylab] = getaxis(a,1); y=double(y);
    [z, zlab] = getaxis(a,3); z=double(z);
    [c, clab] = getaxis(a,0); c=double(c);
    ax = 0; iso = [];
    if ischar(options), 
      iso = str2double(strtok(options));
      if strfind(options, 'axes'), ax=1; end
      if strfind(options, 'auto')
        x=linspace(0,1,size(a,1));
        y=linspace(0,1,size(a,2));
        z=linspace(0,1,size(a,3));
      end
    end
    if (isempty(iso) || ~isfinite(iso))  
      if isnumeric(options) && isscalar(options)
        iso = options;
      else
        iso = (min(c(:))+max(c(:)))/2;
      end
    end
    fv=isosurface(x,y,z,c,iso);
    x3mesh(fv.faces, fv.vertices, 'name', filename, 'subheading', desc, 'rotation', 0, 'axes', ax);
    % Create the supporting X3DOM  files
%    folder = fileparts(filename);
%    load(fullfile(fileparts(which('figure2xhtml')), 'functions', 'x3dom.mat'))
%    folder=fullfile(folder, 'x3dom');
%    if(~isdir(folder)), mkdir(folder); end
%    for i=1:length(data)
%        fid = fopen(fullfile(folder, data(i).filename), 'w','ieee-le');
%        fwrite(fid,data(i).filedata,'uint8');
%        fclose(fid);
%    end
  end
case 'html'
  % create a folder with the HTML style sheet, figures, and document
  [Path, name, ext] = fileparts(filename);
  target = fullfile(Path, name);
  titl = char(a);
  titl(titl=='<')='[';
  titl(titl=='>')=']';
  r = report_generator('iData', target); % directory is filename without extension
  r.open();
  r.section(titl);
  
  % get Model information
  m = []; mp = []; mv = [];
  try
    if isfield(a, 'Model')
      m = get(a, 'Model');
    elseif ~isempty(findfield(a, 'Model'))
      m = get(a, findfield(a, 'Model', 'cache first'));
    end
    
    if isfield(a, 'modelValue')
      mv = get(a, 'modelValue');
    elseif ~isempty(findfield(a, 'modelValue'))
      mv = get(a, findfield(a, 'modelValue', 'cache first'));
    end
    
    if isa(m, 'iFunc')
      % get the parameter values as a struct
      mp = cell2struct(num2cell(m.ParameterValues(:)),strtok(m.Parameters(:)));
    end
  end
  
  % build the HTML content
  r.subsection('Data set')
  % object description
  if ~isempty(a.Title) || ~isempty(title(a))
    data.Title  = strtrim([ a.Title ' ' title(a) ]); end
  if ~isempty(a.Label) || ~isempty(a.DisplayName), 
    data.Label  = strtrim([ a.Label ' ' a.DisplayName ]); end
  data.Source = a.Source;
  data.Date   = [ datestr(a.Date) ', modified ' datestr(a.ModificationDate) ];
  desc = evalc('disp(data)');
  r.add_text([ '<pre> ' desc ' </pre>' ]);
  r.end_section();
  % model description
  if ~isempty(m)
    r.subsection('Model')
    r.add_text(m.Name)
    desc = evalc('disp(mp)');
    r.add_text([ '<pre> ' desc ' </pre>' ]);
    r.end_section();
  end
  
  % plot of the object: special case for 1D which can overlay data and model
  f = [];
  f = figure('Visible','off');
  if ndims(a) == 1 && ~isempty(m) && isa(m, 'iFunc')
    % 1D plot with model
    h=plot(a,'bo',m,'r-','tight');
  elseif ndims(a) == 1
    % simple plot. Model not available.
    h=plot(a,'bo','tight');
  elseif ndims(a) > 1
    if ~isempty(m) && isa(m, 'iFunc')
      h=subplot([a m], [1 2], 'tight');
    else
      h=plot(a, 'tight');
    end
  end
  if ishandle(f)
    set(f, 'Name', [ 'iFit_DataSet_' a.Tag ]);
    r.add_figure(f, char(a), 'centered');
  end
  close(f);
  % export object into a number of usable formats
  builtin('save', fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.mat' ]), 'a');
  save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.pdf' ]), 'pdf');
  save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.fig' ]), 'fig');
  save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.h5' ]), 'mantid');
  if prod(size(a)) < 1e5
    save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.dat' ]), 'dat data');
    save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.svg' ]), 'svg');
  end
  t = 'Exported to: [ ';
  for ext={'mat','png','dat','svg','pdf','fig','h5'}
    f = [ 'iFit_DataSet_' a.Tag '.' ext{1} ];
    if ~isempty(dir(fullfile(target, 'img', f)))
      t = [ t '<a href="img/' f '"> ' ext{1} ' </a>' ];
    end
  end
  t = [ t ' ]' ];
  r.add_text(t);
  
  % add more information about the Data set
  r.subsection('More about the Data set...')
  desc = evalc('disp(a,''data'',''flat'')');
  desc(desc=='<')='[';
  desc(desc=='>')=']';
  r.add_text([ '<br><pre> ' desc ' </pre>' ]);
  
  % display a 'footer' below the object description
  r.add_text('<hr>');
  r.add_text([ '<b>' datestr(now) '</b> - ' version(iData) '<br>' ])
  
  r.add_text([ '<a href="http://ifit.mccode.org">Powered by iFit ' ...
    '<img src="http://ifit.mccode.org/images/iFit-logo.png" width=35 height=32></a> ' ...
    '<a href="http://www.ill.eu">(c) ILL ' ...
    '<img title="ILL, Grenoble, France www.ill.eu" src="http://ifit.mccode.org/images/ILL-web-jpeg.jpg" alt="ILL, Grenoble, France www.ill.eu" style="width: 21px; height: 20px;"></a><hr>' ]);
  r.end_section();
  r.close();
case {'stl','stla','stlb','off','ply'} % STL ascii, binary, PLY, OFF
  if ndims(a) == 1    iData_private_warning(mfilename,[ 'Object ' inputname(1) ' ' a.Tag ' 1D does not seem to be exportatble as a ' format ' file. Ignoring.' ]);
    return
  else
    if ndims(a) == 2
      [x, xlab] = getaxis(a,2); x=double(x);
      [y, ylab] = getaxis(a,1); y=double(y);
      [z, zlab] = getaxis(a,0); z=double(z);
    elseif ndims(a) >= 3
      [x, xlab] = getaxis(a,2); x=double(x);
      [y, ylab] = getaxis(a,1); y=double(y);
      [z, zlab] = getaxis(a,3); z=double(z);
    else
      disp([ mfilename ': Export into ' format ' is only possible for 2D+ objects, not for ' num2str(ndims(a)) 'D. Ignoring.' ])
      return
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
    if isvector(x) && isvector(y)
      [x,y] = meshgrid(x,y);
    end
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
case {'yaml','yml'}
  YAML.write( filename, struct(a) ); % YAML object is in iFit/Objects
case 'json'
  mat2json(struct(a), filename );    % in private
case {'xml'}
  try
      struct2xml(struct(a), filename);   % in private
  catch
      filename = [];
  end
otherwise
  iData_private_warning(mfilename,[ 'Export of object ' inputname(1) ' ' a.Tag ' into format ' format ' is not supported. Ignoring.' ]);
  filename = [];
end

% end of iData/saveas



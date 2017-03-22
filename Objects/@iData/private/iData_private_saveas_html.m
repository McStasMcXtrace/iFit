function filename = iData_private_saveas_html(a, filename, format)
  % save an iData into an HTML document
  % if the document already exists, the new content is appended
  %
  % 'data' in format triggers minimalistic mode for printing.
  if isempty(a), filename=[]; return; end
  if isempty(dir(filename))
    mode = 'w+';
  else 
    mode = 'a+';
  end
  if nargin < 3, format = ''; else format=lower(format); end

  [Path, name, ext] = fileparts(filename);
  target = Path;
  titl = char(a);
  titl(titl=='<')='[';
  titl(titl=='>')=']';
  titl(titl=='\')='';
  % Open and write the HTML header
  fid = fopen(filename, mode);  % create or append to file
  if fid == -1, filename = []; return; end
  if ~isdir(fullfile(target,'img')), mkdir(fullfile(target,'img')); end
  
  % The Header *****************************************************************
  if strcmp(mode, 'w+')
    fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n');
    fprintf(fid, '<html>\n<head>\n<title>%s</title>\n</head>\n', ...
        titl);
    fprintf(fid, '<body>\n');
  end
  fprintf(fid, '<div style="text-align: center;"><h1>%s</h1></div>\n', titl);
  
  % get any Model information
  m = []; mp = []; mv = [];
  try
    if isfield(a, 'Model')
      m = get(a, 'Model');
    elseif ~isempty(findfield(a, 'Model'))
      m = get(a, findfield(a, 'Model', 'cache first'));
    end
    
    if isfield(a, 'ModelValue')
      mv = get(a, 'ModelValue');
    elseif ~isempty(findfield(a, 'ModelValue'))
      mv = get(a, findfield(a, 'ModelValue', 'cache first'));
    end
    
    if isfield(a, 'ModelParameters')
      mp = get(a, 'ModelParameters');
      if isa(m, 'iFunc') && ~isempty(m)
        names = m.Parameters(:);
      end
      if numel(names) == numel(mp)
        mp = cell2struct(num2cell(mp(:)),strtok(names(:)));
      end
    elseif isfield(a, 'Parameters')
      mp = get(a, 'Parameters');
    elseif isa(m, 'iFunc') && ~isempty(m)
      % get the parameter values as a struct
      mp = cell2struct(num2cell(m.ParameterValues(:)),strtok(m.Parameters(:)));
    end
  end
  
  % build the HTML content
  
  % Data set information (short) ***********************************************
  fprintf(fid,'<h2>Data set</h2>\n');
  % object description
  data = [];
  if ~isempty(a.Title) || ~isempty(title(a))
    data.Title  = strtrim([ a.Title ' ' title(a) ]); end
  if ~isempty(a.Label) || ~isempty(a.DisplayName), 
    data.Label  = strtrim([ a.Label ' ' a.DisplayName ]); end
  data.Source = a.Source;
  data.Date = datestr(now);
  if isnumeric(a.Date)
    data.Date   = datestr(a.Date);
  elseif ischar(a.Date)
    data.Date = a.Date;
  end
  if isnumeric(a.ModificationDate)
    data.Date   = [ data.Date ', modified ' datestr(a.ModificationDate) ];
  elseif ischar(a.ModificationDate)
    data.Date   = [ data.Date ', modified ' a.ModificationDate ];
  end
  desc = evalc('disp(data)');
  fprintf(fid,[ '<pre> ' desc ' </pre>\n' ]);

  % Model information **********************************************************
  % model description when exists
  if ~isempty(mp)
    desc = evalc('disp(mp)');
  else desc = ''; end
  
  if ~isempty(m)
    fprintf(fid,'<h2>Model: %s</h2>\n', m.Name);
    fprintf(fid, '<p>\n');
    fprintf(fid, m.Description);
    fprintf(fid, '</p>\n');
  end
  if ~isempty(mp)
    fprintf(fid, 'Model parameters:<br>\n');
    fprintf(fid, [ '<pre>' desc '</pre>\n' ]);
  end
  
  % Data set (and Model) plot **************************************************
  % plot of the object: special case for 1D which can overlay data and model
  f = figure('Visible','off', 'Name', [ 'iFit_DataSet_' a.Tag ]);
  if ndims(a) == 1
    % simple plot and model
    h=plot(a,'bo','tight');
  elseif ndims(a) > 1
    if ndims(a) == 2 || isempty(mv)
      % overlay data and model
      h = plot(a, 'tight view3');
    else
      % side by side
      h = subplot([a mv], [1 2], 'tight');
    end
  end
  % add text with parameters onto plot
  if ~isempty(desc)
    if ~isempty(m)
      desc = sprintf('%s\n%s', m.Name, desc);
    end
    h = text(0,0, desc, 'Unit','normalized','Interpreter','none', ...
      'BackgroundColor',[0.9 0.9 0.9],'FontName','FixedWidth');
  end
  
  % Export data ****************************************************************
  
  % create output from the figure: png pdf fig
  basename     = fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag ]);
  basename_img = fullfile('img', [ 'iFit_DataSet_' a.Tag ]);
  if ndims(a) == 2
    saveas(a, basename, 'png data');  % just the image for 2D data sets, else getframe
  elseif ndims(a) == 3
    saveas(a, basename, 'png', 'tight view3');
  else
    saveas(a, basename, 'png', 'tight');
  end
  
  % only export when not in 'data' mode (flat HTML with minimal content)
  if isempty(strfind(format, 'data')) && isempty(strfind(format, 'flat'))
    saveas(f, basename, 'fig');
    saveas(f, basename, 'pdf');
  end
  close(f);
  
  % export object into a number of usable formats
  % 1D: mat dat/data hdf5/mantid png json      xml      yaml      nc fig pdf svg
  % 2D: mat dat/data hdf5/mantid png json      xml      yaml      nc fig pdf svg  fits
  % 3D: mat dat/data hdf5/mantid png json/data xml/data yaml/data nc fig pdf vtk mrc nc   
  % nD: mat dat/data hdf5/mantid png json/data xml/data yaml/data nc
  export       = {'mat','dat data','hdf mantid', 'json','xml','yaml','nc','tiff' };
  export_label = { ...
  'Matlab binary file. Open with Matlab or <a href="http://ifit.mccode.org">iFit</a>.', ...
  'Flat text file which contains axes and the data set. You will have to reshape the matrix after reading the contents. View with any text editor.', ...
  '<a href="http://www.hdfgroup.org/">NeXus/HDF5</a> data file, to be opened with e.g. <a href="http://www.mantidproject.org/Main_Page">Mantid</a>, <a href="http://www.hdfgroup.org/hdf-java-html/hdfview">hdfview</a> or <a href="http://ifit.mccode.org">iFit</a>.', ...
  '<a href="http://en.wikipedia.org/wiki/JSON">JavaScript Object Notation</a>, to be opened with e.g. JSONView Chrome/Firefox plugin and text editors.', ...
  '<a href="http://www.w3.org/XML/">Extensible Markup Language</a> file, to be opened with e.g. Chrome/Firefox and text editors.', ...
  '<a href="http://en.wikipedia.org/wiki/YAML">YAML</a> interchange format, to be viewed with e.g. text editors.', ...
  '<a href="http://www.unidata.ucar.edu/software/netcdf/">NetCDF</a> binary file, to be viewed with <a href="http://meteora.ucsd.edu/~pierce/ncview_home_page.html">ncview</a> and <a href="http://www.hdfgroup.org/hdf-java-html/hdfview">hdfview</a>.', ...
  '<a href="https://en.wikipedia.org/wiki/TIFF">TIFF</a> image file, to be viewed with e.g. <a href="http://rsb.info.nih.gov/ij/">ImageJ</a>, <a href="http://www.gimp.org/">GIMP.</a>.'};
  
  % add 'data' keyword when the object memory size is larger than 10 Mb
  w = whos('a');
  flag_data  = (w.bytes > 10e6);
  
  if ~flag_data
    export = [ export 'svg' ];  % when not too big
    export_label = [ export_label '<a href="https://fr.wikipedia.org/wiki/Scalable_Vector_Graphics">Scalable Vector Graphics</a> image, to be viewed with Chrome/Firefox, <a href="http://inkscape.org/">Inkscape</a>, <a href="http://www.gimp.org/>GIMP.</a>, <a href="http://projects.gnome.org/evince/">Evince</a>.' ];
  end
  if ndims(a) == 3
    export = [ export 'vtk' 'mrc' ];
    export_label = [ export_label ...
      'Visualization Toolkit (VTK) file which can be viewed with <a href="http://www.paraview.org/">ParaView</a>, <a href="http://code.enthought.com/projects/mayavi/">Mayavi2</a>, <a href="https://wci.llnl.gov/simulation/computer-codes/visit/executables">VisIt</a>, <a href="https://www.slicer.org/">Slicer4</a>.', ...
      'MRC Electron density map, to be visualized with <a href="http://www.pymol.org/">PyMol</a>, <a href="http://www.ks.uiuc.edu/Research/vmd/">VMD</a>, <a href="http://www.cgl.ucsf.edu/chimera/">Chimera</a>, <a href="http://www.yasara.org/">Yasara</a>, <a href="http://mem.ibs.fr/VEDA/">VEDA</a>.' ];
  end
  if ~flag_data && any(ndims(a) == [2 3])
    export = [ export 'xhtml' 'x3d' ];
    export_label = [ export_label ...
      'Extensible Web page with embeded viewer (X3DOM), to be viewed with Chrome/Firefox.', ...
      'X3D Geometry Scene for <a href="http://castle-engine.sourceforge.net/view3dscene.php">view3dscene</a>, <a href="http://www.instantreality.org/">InstantPlayer</a>, <a href="http://freewrl.sourceforge.net/">FreeWRL</a>' ];
  end
  if ndims(a) == 2
    export = [ export 'fits' ];
    export_label = [ export_label ...
      'Flexible Image Transport System (<a href="https://fits.gsfc.nasa.gov/fits_home.html">FITS</a>) image, which can be viewed with e.g. <a href="http://rsb.info.nih.gov/ij/">ImageJ</a>, <a href="http://www.gimp.org/">GIMP.</a>.' ];
  end
  if isempty(strfind(format, 'data')) && isempty(strfind(format, 'flat'))
    for index=1:numel(export)
      f = export{index};
      if flag_data && isempty(strfind(f, 'data'))
        f = [ f ' data' ];
      end
      switch export{index}
      case 'mat'
        builtin('save', basename, 'a');
      otherwise
        save(a, basename, f);
      end
    end
  end
  export = [ export 'png' 'fig' 'pdf' ];
  export_label = [ export_label, ...
    'PNG image for <a href="http://www.gimp.org/">GIMP</a> or <a href="http://projects.gnome.org/evince/">Evince</a>', ...
    'Matlab figure to be opened with Matlab or <a href="http://ifit.mccode.org">iFit</a>. Use <i>set(gcf,''visible'',''on'')</i> after loading.', ...
    'Portable Document File to be viewed with <a href="http://get.adobe.com/fr/reader/">Acrobat Reader</a> or <a href="http://projects.gnome.org/evince/">Evince</a>.' ];
  
  % add image and links to exported files
  if ~isempty(dir([ basename '.png' ]))
    fprintf(fid, '<div style="text-align: center;"><a href="%s"><img src="%s" align="middle"></a><br>\n<i>Data: %s</i><br></div>\n', ...
      [ basename_img '.png' ], ...
      [ basename_img '.png' ], titl);
    if ndims(a) == 2
      fprintf(fid, '(try the <a href="%s">TIFF file</a> in case the axes are not shown)<br>\n', [ basename_img '.tiff' ]);
    end
    if ~isempty(m)
      fprintf(fid, '<div style="text-align: center;"><i>Model: %s</i><br></div>\n', m.Name);
    end
  end
  
  % display list of available formats, as well as suggested software to use
 fprintf(fid, '<p>Exported to: <br><ul>\n');
  for index=1:numel(export)
    if ~isempty(dir([ basename '.' strtok(export{index}) ]))
      fprintf(fid, [ '<li><b><a href="' basename_img '.' strtok(export{index}) '">' export{index} '</a></b>: ' ...
        export_label{index} '</li>\n' ]);
    elseif isempty(strfind(format, 'data'))
      disp([ mfilename ': skipping ' basename '.' strtok(export{index}) ])
    end
  end
  fprintf(fid, '</ul></p>\n');
  % special case for XHTML (displayed in a frame)
  if ~isempty(dir([ basename '.' 'xhtml' ]))
    fprintf(fid, [ '<iframe src="%s" align="middle" width="1024" height="1024"></iframe><br>' ...
    '(<a href="%s" target=_blank>open in external window</a>)<br>\n' ], ...
    [ basename_img '.xhtml' ], [ basename_img '.xhtml' ]);
  end
  
  % Data set information (details) *********************************************
  fprintf(fid,'<h2>More about the Data set...</h2\n');
  desc = evalc('disp(a,''data'',''flat'')');
  desc(desc=='<')='[';
  desc(desc=='>')=']';
  fprintf(fid,[ '<br><pre> ' desc ' </pre>\n' ]);
  fprintf(fid,'<hr>\n');
  
  % The Footer *****************************************************************
  if strcmp(mode, 'w+')
    % display a 'footer' below the object description
    
    fprintf(fid,[ '<b>' datestr(now) '</b> - ' version(iData) '<br>\n' ]);
    
    fprintf(fid,[ '<a href="http://ifit.mccode.org">Powered by iFit ' ...
      '<img src="http://ifit.mccode.org/images/iFit-logo.png" width=35 height=32></a> \n' ...
      '<a href="http://www.ill.eu">(c) ILL ' ...
      '<img title="ILL, Grenoble, France www.ill.eu" src="http://ifit.mccode.org/images/ILL-web-jpeg.jpg" alt="ILL, Grenoble, France www.ill.eu" style="width: 33px; height: 32px;"></a><hr>\n' ]);
  end
  
  fprintf(fid,'<p><!-- pagebreak --></p>\n'); % force page break in case we append new stuff
  fclose(fid);


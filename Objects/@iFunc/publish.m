function filename = publish(a, filename, section, message)
% b = publish(s) : export the Model object as an HTML document.
%
%   @iFunc/publish function to export an object into a readable document.
%
%   publish(a)        export to a temporary directory and open the web browser
%   publish(a, file)  export to given file name
%   publish(a, dir)   export to given directory
%   publish(a, file, sections) export only the given sections as char or cell 
%     default sections are {'header','parameters','model','expression','footer'};
%     optional section:     'download','message'
%   publish(a, file, section, message)
%     write a specified section with additional label/text in 'message'.
%   publish(a, file, 'force')
%     overwrites any existing document (default is to append)
%
% When the generated report file already exists, the new document is appended.
% The resulting document file name is stored as:
%   a.UserData.options.report
%
% input:  a: object or array (iFunc)
%         filename: a file name or directory where to export. When not given
%           the exportation takes place in the local directory.
%         section: a list of sections to generate. The default is
%           {'header','plot','expression','footer'}. 
%            'download' and 'message' can be added.
%         message: an additional label for the section.
% output: b: generated filename
% ex:     publish(gauss)
%
% Version: $Date$
% See also iFunc, iFunc/save

if nargin < 2, filename = ''; end
if nargin < 3, section  = []; end
if nargin < 4, message  = []; end
if strcmp(filename,'force') && isempty(section)
  filename=''; section='force';
end

% handle filename and target directory
UserData = a.UserData; % this is where we search for information
if isfield(UserData, 'options'), options = UserData.options; else options = []; end

% determine where to store the results =======================================
% output: is there a location already ?
if isfield(options, 'report')
  filename = options.report;
else
  filename = findfield(a, 'target');
  if isempty(filename)    filename = findfield(a, 'filename'); end
  if isempty(filename)    filename = findfield(a, 'report'); end
  if isempty(filename)    filename = findfield(a, 'dir'); end
  if iscell(filename)     filename = filename{1}; end
  if ~isempty(filename)   filename = get(a,filename); end
end
% output: add extension when not set
if isdir(char(filename))
  p = filename; f = []; e = [];
else
  [p,f,e] = fileparts(char(filename));
end
if isempty(p)  p = tempname; end
if isempty(f), f = [ 'iFit_Model_' a.Tag ]; end
if isempty(e), e = '.html'; end
filename       = fullfile(p, [f e ]);  % store location and file name
options.report = filename;
options.target = p;

% overwrite when in 'force' mode
if any(strcmp(section, 'force'))
  if ~isempty(dir(filename)), delete(filename); end
  if iscell(section), section(strcmp(section,'force')) = [];
  else section = []; end
end

% create directories if needed
if ~isdir(options.target), mkdir(options.target); end

% now write the documents into bits
a = publish_write(a, filename, section, message);

% finalize: store back all into current object
UserData.options = options;
a.UserData       = UserData;
if ~isempty(inputname(1))
  assignin('caller',inputname(1),a); % update in original object
end

% last, open web browser when no output requested ------------------------------
if nargout == 0
  webbrowser(filename, 'system');
end





% ------------------------------------------------------------------------------
function a = publish_write(a, filename, section, message)
  % publish_write: actually write sections

  if nargin < 3, section = []; end
  if isempty(section)
    section = {'header','parameters','plot','model','expression','footer'};
  end
  
  if iscell(section)
    for index=1:numel(section)
      publish_write(a, filename, section{index}, message);
      if numel(section) > 1, message=[]; end  % message only displayed once
    end
    return
  end

  if isempty(dir(filename))
    mode = 'w+';
  else 
    mode = 'a+';
  end
  target = fileparts(filename);
  titl   = a.Name;
  titl(titl=='<')='[';
  titl(titl=='>')=']';
  
  % Open and write the HTML header
  fid = fopen(filename, mode);  % create or append to file
  if fid == -1, filename = []; return; end
  if ~isdir(fullfile(target,'img')), mkdir(fullfile(target,'img')); end
  
  basename     = fullfile(target, 'img', [ 'iFit_Model_' a.Tag ]);
  basename_img = fullfile('img',         [ 'iFit_Model_' a.Tag ]);
  
  Section = lower(strtok(section)); Section(1) = upper(Section(1));
  
  switch strtok(lower(section))
  
  
  case {'header','force'}
  
    % The Header ***************************************************************
    % Model description
    if strcmp(mode, 'w+')
      fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n');
      fprintf(fid, '<html>\n<head>\n<title>%s</title>\n</head>\n', ...
          titl);
      fprintf(fid, '<body>\n');
    end
    fprintf(fid, '<h1>Model: %s</h1>\n<p><pre>%s</pre></p>\n', titl, a.Description);
    fprintf(fid, '<b>Date</b>: %s<br>\n<b>Tag:</b>: %s<br>\n', datestr(now), a.Tag);
    fprintf(fid, '<b>Class</b>: %s</br>\n', class(a));
    fprintf(fid, '<b>Stored in</b>: <a href="%s">%s</a><br>\n', target, target);
    
  case 'parameters'
  
    % Model Parameters
    fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', section, Section);
    
    fprintf(fid, 'Model "%s" [%s] parameters:<br><br><table>\n', a.Name, a.Tag); % displayed as a list
    for index=1:numel(a.Parameters)
      [par, rem] = strtok(a.Parameters{index}); % get parameter name and possibly any other comment
      if ~isempty(a.ParameterValues)
        fprintf(fid, '<tr><td><b>%s</b> %s</td><td>%g</td></tr>\n', par, rem, a.ParameterValues(index));
      else
        fprintf(fid, '<tr><td><b>%s</b> %s</td></tr>\n', par, rem);
      end
    end
    fprintf(fid,'</table>\n');
    
  
  case 'plot'

    % Model plot ***************************************************************
    
    % get parameter values
    if ~isempty(a.ParameterValues)
      mp = cell2struct(num2cell(a.ParameterValues(:)),strtok(a.Parameters(:)));
      desc = evalc('disp(mp)');
    else
      mp = a.Parameters;
      desc = ''; 
    end
    
    f = figure('Visible','off', 'Name', [ 'iFit_Model_' a.Tag ]);
    h = plot(a); axis tight;
    % add text with parameters onto plot
    if ~isempty(desc)
      desc = sprintf('%s\n%s', titl, desc);
      h = text(0,0, desc, 'Unit','normalized','Interpreter','none', ...
        'BackgroundColor',[0.9 0.9 0.9],'FontName','FixedWidth');
    end

    % create output from the figure: png pdf fig
    set(f,'Renderer','zbuffer');
    saveas(f, basename, 'fig');
    saveas(f, basename, 'png');
    saveas(f, basename, 'pdf');
    saveas(f, basename, 'epsc');
    plot2svg([ basename '.svg' ], f, 'jpg');
    close(f);
    
  case 'model'
    
    % export object into a number of usable formats
    fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', section, Section);
    export       = {'mat','dat','json','xml','yaml'};
    
    for index=1:numel(export)
      f = export{index};
      switch f
      case 'mat'
        builtin('save', basename, 'a');
      otherwise
        save(a, basename, f);
      end
    end
    
    % display list of available Model formats, as well as suggested software to use
    fprintf(fid, '<h4>Model "%s" [%s] Exported to: </h4><table>\n', titl, a.Tag);
    publish_table(fid, target, 'img', [ 'iFit_Model_' a.Tag ]);
  
  case 'expression'
  
    % Model expression (details)
    fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', section, Section);
    
    t = cellstr(a);
    fprintf(fid, '<pre>\n');
    fprintf(fid, '%s\n', t{:});
    fprintf(fid, '</pre>\n');
    
  case {'status','message'}
  
    % catenate a single message (HTML formatting)
    if ~isempty(message)
      if iscell(message)
        fprintf(fid, '%s\n', message{:});
      elseif isstruct(message)
        % write a Table
      elseif ischar(message)
        fprintf(fid, '%s\n', message);
      end
      message = [];
    end
  
  case 'download'
  
    % create a ZIP of the 'target' directory and display its link for download
    fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', section, Section);
    
    fprintf(fid, 'You can download the whole content of this report (data, plots, ...) from<br>\n');
    [~,short_path] = fileparts(target);
    fprintf(fid, '<ul><li><a href="%s">%s</a></li></ul>', fullfile('..',[ short_path '.zip' ]), [ short_path '.zip' ]);
    fprintf(fid, 'You should then open the "html" file therein with a browser (Firefox, Chrome...)<br><br>\n');
    
    % create a simple README file
    freadme = fopen(fullfile(target,'README.txt'),'w');
    fprintf(freadme, 'Model: %s\n%s\n\n', titl, a.Description);
    fprintf(freadme, [ 'Open the .html file in this directory.\n', ...
      'It contains all you need. ', ...
      'The Model (%s) is also stored in the .mat Matlab binary file.\n\n' ], class(a));
    fprintf(freadme, [ datestr(now) ' - ' version(iData) '\n' ]);
    fclose(freadme);
    
    % create ZIP of document
    zip(target, target);
    disp([ mfilename ': compressed as ' target '.zip' ]);
  
  case {'footer'}
  
    % The Footer ***************************************************************
    % display a 'footer' below the object description
    fprintf(fid,[ '<hr><br><b>' datestr(now) '</b> - ' version(iData) '<br>\n' ]);
    
    fprintf(fid,[ '<a href="http://ifit.mccode.org">Powered by iFit ' ...
      '<img src="http://ifit.mccode.org/images/iFit-logo.png" width=35 height=32></a> \n' ...
      '<a href="http://www.ill.eu">(c) ILL ' ...
      '<img title="ILL, Grenoble, France www.ill.eu" src="http://ifit.mccode.org/images/ILL-web-jpeg.jpg" alt="ILL, Grenoble, France www.ill.eu" style="width: 33px; height: 32px;"></a><hr>\n' ]);
  
    fprintf(fid,'<p><!-- pagebreak --></p>\n'); % force page break in case we append new stuff
    
    if isdeployed || ~usejava('jvm') || ~usejava('desktop')
      disp([ mfilename ': HTML report created as ' filename ]);
    else
      disp([ mfilename ': HTML report created as <a href="' filename '">' filename '</a>' ]);
    end
  case 'table'
    if ischar(message)
      publish_table(fid, target, '', message);
      publish_table(fid, target, 'img', message);
      message = [];
    end
  end % switch
  
  if ~isempty(message), publish(self, filename, 'message', message); end
  fclose(fid);


% ------------------------------------------------------------------------------
function publish_table(fid, target, img, name)
% write a Table which displays available exported 'configuration' files
  
  if fid == -1, return; end
  if isempty(name), name = 'configuration'; end
  if isempty(img),  img = ''; end
  % append links to configuration files and images
  flag_table = 0;
  for index={ ...
    '.html You can rotate the object (left mouse button), zoom (right mouse button), and pan (middle mouse button)', ...
  '.png', ...
  '.cif Crystallographic Information File (<a href="https://en.wikipedia.org/wiki/Crystallographic_Information_File">CIF</a>) which you can view with <a href="http://jmol.sourceforge.net/">JMol</a>, <a href="http://www.ks.uiuc.edu/Research/vmd/">VMD</a>, <a href="http://www.cgl.ucsf.edu/chimera/">Chimera</a>, <a href="http://rasmol.org/">RasMol</a>, ...', ...
  '.etsf <a href="http://www.etsf.eu/fileformats">European Theoretical Spectroscopy Facility</a> file format for DFT codes', ...
  '.pdb Protein Data Bank (<a href="http://www.rcsb.org/pdb/">PDB</a>) which you can view with <a href="http://jmol.sourceforge.net/">JMol</a>, <a href="http://www.ks.uiuc.edu/Research/vmd/>"VMD</a>, <a href="http://www.cgl.ucsf.edu/chimera/">Chimera</a>, <a href="http://rasmol.org/">RasMol</a>, ...', ...
  '_POSCAR POSCAR geometry file for <a href="https://www.vasp.at/">VASP</a>', ...
  '_SHELX.res <a href="http://shelx.uni-ac.gwdg.de/SHELX/">ShelX</a> file format', ...
  '.eps Encapsulated postscript', ...
  '.pov Scene for <a href="http://www.povray.org/">Pov-Ray</a>', ...
  '.x3d Geometry Scene for <a href="http://castle-engine.sourceforge.net/view3dscene.php">view3dscene</a>, <a href="http://www.instantreality.org/">InstantPlayer</a>, <a href="http://freewrl.sourceforge.net/">FreeWRL</a>' ...
  '.json <a href="http://en.wikipedia.org/wiki/JSON">JavaScript Object Notation</a>, to be opened with e.g. JSONView Chrome/Firefox plugin and text editors.' ...
  '.xml <a href="http://www.w3.org/XML/">Extensible Markup Language</a> file, to be opened with e.g. Chrome/Firefox and text editors.' ...
  '.yaml <a href="http://en.wikipedia.org/wiki/YAML">YAML</a> interchange format, to be viewed with e.g. text editors.' ...
 [ '.mat Matlab binary file. Load the Model/Data set under <a href="http://ifit.mccode.org">Matlab/iFit</a> with (this also works with the <a href="http://ifit.mccode.org/Install.html">standalone version of iFit</a> which does <b>not</b> require any Matlab license and installs on most systems): load(''<a href="' name '.mat">' name '.mat</a>'') <i></i></li></ul>' ] ...
  '.m Matlab script/function' ...
  '.dat Flat text file. You may have to reshape the data set. View with any text editor (gedit).' ...
  '.h5 a NeXus/HDF5 data file to be opened with e.g. <a href="http://www.mantidproject.org/Main_Page">Mantid</a>, <a href="http://www.hdfgroup.org/hdf-java-html/hdfview">hdfview</a>, <a href="http://vitables.org/">ViTables</a> or <a href="http://ifit.mccode.org">iFit</a>' ...
  '.fig a Matlab figure for Matlab or <a href="http://ifit.mccode.org">iFit</a>. Use </i>figure(gcf)</i> after loading' ...
  '.pdf an Adobe PDF, to be viewed with <a href="http://get.adobe.com/fr/reader/">Acrobat Reader</a> or <a href="http://projects.gnome.org/evince/">Evince</a>' ...
  '.svg a <a href="https://fr.wikipedia.org/wiki/Scalable_Vector_Graphics">Scalable Vector Graphics</a> image, to be viewed with Chrome/Firefox, <a href="http://inkscape.org/">Inkscape</a>, <a href="http://www.gimp.org/>GIMP.</a>, <a href="http://projects.gnome.org/evince/">Evince</a>' ...
  '.vtk Visualization Toolkit (VTK) file which can be viewed with <a href="http://www.paraview.org/">ParaView</a>, <a href="http://code.enthought.com/projects/mayavi/">Mayavi2</a>, <a href="https://wci.llnl.gov/simulation/computer-codes/visit/executables">VisIt</a>, <a href="https://www.slicer.org/">Slicer4</a>' ...
  '.mrc MRC Electron density map, to be visualized with <a href="http://www.pymol.org/">PyMol</a>, <a href="http://www.ks.uiuc.edu/Research/vmd/">VMD</a>, <a href="http://www.cgl.ucsf.edu/chimera/">Chimera</a>, <a href="http://www.yasara.org/">Yasara</a>, <a href="http://mem.ibs.fr/VEDA/">VEDA</a>' ...
  '.xhtml Extensible Web page. You can rotate the object (left mouse button), zoom (right mouse button), and pan (middle mouse button).', ...
  '.fits Flexible Image Transport System (<a href="https://fits.gsfc.nasa.gov/fits_home.html">FITS</a>) image, which can be viewed with e.g. <a href="http://rsb.info.nih.gov/ij/">ImageJ</a>, <a href="http://www.gimp.org/">GIMP.</a>.', ...
  '.tiff <a href="https://en.wikipedia.org/wiki/TIFF">TIFF</a> image file, to be viewed with e.g. <a href="http://rsb.info.nih.gov/ij/">ImageJ</a>, <a href="http://www.gimp.org/">GIMP.</a>.' ...
    }
    
    [index1, index2] = strtok([ name index{1} ]);
    if ~isempty(dir(fullfile(target,img,index1)))
      % the file exists
      if strcmp(index1, [ name '.png' ]) && ~isempty(dir(fullfile(target,img,[ name '.tiff' ])))
        fprintf(fid, [ '<div style="text-align:center">' ...
          '<a href="%s"><img src="%s" align="middle" title="%s" height="480" width="640">' ...
          '</a><br>(try the <a href="%s">TIFF file</a>' ...
          ' in case the axes are not shown)</div><br>\n' ], ...
          fullfile(img,index1), fullfile(img,index1), fullfile(img,index1), [ name '.tiff' ]);
      elseif strcmp(index1, [ name '.png' ])  % only the PNG, but not the TIFF
        fprintf(fid, [ '<div style="text-align:center">' ...
          '<a href="%s"><img src="%s" align="middle" title="%s" height="480" width="640">' ...
          '</a></div><br>\n' ], fullfile(img,index1), fullfile(img,index1), fullfile(img,index1));
      elseif strcmp(index1, [ name '.html' ])
        % embed a frame
        fprintf(fid, [ '<div style="text-align:center">' ...
          '<iframe src="%s" align="middle" width="480" height="480">' ...
          '</iframe><br>%s<br>' ...
          '(<a href="%s" target=_blank>open in external window</a>)<br></div><br>\n' ], ...
          fullfile(img,index1), index2, fullfile(img,index1));
      elseif strcmp(index1, [ name '.xhtml' ])
        % embed an x3dom frame
        fprintf(fid, [ '<div style="text-align:center"><iframe src="%s" align="middle" width="700" height="850"></iframe><br>\n' ...
        '%s<br>(<a href="%s" target=_blank>open in external window</a>)<br>\n' ...
        'You can rotate the model (left mouse button), zoom (right mouse button), and pan (middle mouse button).<br>\n' ], fullfile(img,index1), index2, fullfile(img,index1));
      else
        if flag_table == 0
          flag_table = 1;
          fprintf(fid, '<table>\n');
        end
        fprintf(fid, '  <tr><td><a href="%s">%s</a></td><td>%s</td></tr>\n', fullfile(img,index1), index1, index2);
      end
    end
  end
  if flag_table
    fprintf(fid, '</table><br>\n');
  end

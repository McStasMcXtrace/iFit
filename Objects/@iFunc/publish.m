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
  if iscell(filename)     filename = filename{1}; end
  if ~isempty(filename)   filename = get(a,filename); end
end
% output: add extension when not set
[p,f,e] = fileparts(char(filename));
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
publish_write(a, filename, section, message);

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
function publish_write(a, filename, section, message)
  % publish_write: actually write sections

  if nargin < 3, section = []; end
  if isempty(section)
    section = {'header','parameters','model','expression','footer'};
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
  
  % get parameter values
  if ~isempty(a.ParameterValues)
    mp = cell2struct(num2cell(a.ParameterValues(:)),strtok(a.Parameters(:)));
    desc = evalc('disp(mp)');
  else
    mp = a.Parameters;
    desc = ''; 
  end
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
    fprintf(fid, '<h1>Model: %s</h1>\n<p>%s</p>\n', titl, a.Description);
    fprintf(fid, '<b>Date</b>: %s<br>\n<b>Computer</b>: %s<br>\n', datestr(now), computer);
    fprintf(fid, '<b>Class</b>: %s</br>\n', class(a));
    fprintf(fid, '<b>Stored in</b>: <a href="%s">%s</a><br>\n', target, target);
    
  case 'parameters'
  
    % Model Parameters
    fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', Section, Section);
    
    fprintf(fid, 'Model parameters<br>\n'); % displayed as a list
    fprintf(fid,'<ul>\n');
    for index=1:numel(a.Parameters)
      if ~isempty(a.ParameterValues)
        fprintf(fid, '<li>%s = %g</li>\n', a.Parameters{index}, a.ParameterValues(index));
      else
        fprintf(fid, '<li>%s</li>\n', a.Parameters{index});
      end
    end
    fprintf(fid,'</ul>\n');
    
  
  case 'model'

    % Model plot ***************************************************************
    % exported plots/files
    fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', Section, Section);
    
    f = figure('Visible','off', 'Name', [ 'iFit_Model_' a.Tag ]);
    h = plot(a); axis tight;
    % add text with parameters onto plot
    if ~isempty(desc)
      desc = sprintf('%s\n%s', a.Name, desc);
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

    % export object into a number of usable formats
    export       = {'mat','dat',' hdf4', 'json','xml','yaml'};
    export_label = { ...
    'Matlab binary file. Open with Matlab or <a href="http://ifit.mccode.org">iFit</a>.', ...
    'Flat text file which contains the Model expression', ...
    '<a href="http://www.hdfgroup.org/">HDF4</a> image, to be opened with e.g. <a href="http://www.hdfgroup.org/hdf-java-html/hdfview">hdfview</a> or <a href="http://ifit.mccode.org">iFit</a>.', ...
    '<a href="http://en.wikipedia.org/wiki/JSON">JavaScript Object Notation</a>, to be opened with e.g. JSONView Chrome/Firefox plugin and text editors.', ...
    '<a href="http://www.w3.org/XML/">Extensible Markup Language</a> file, to be opened with e.g. Chrome/Firefox and text editors.', ...
    '<a href="http://en.wikipedia.org/wiki/YAML">YAML</a> interchange format, to be viewed with e.g. text editors.' };
    
    for index=1:numel(export)
      f = export{index};
      switch f
      case 'mat'
        builtin('save', basename, 'a');
      otherwise
        save(a, basename, f);
      end
    end
    export = [ export 'png' 'fig' 'pdf' 'svg' 'eps' ];
    export_label = [ export_label, ...
      'PNG image for <a href="http://www.gimp.org/">GIMP</a> or <a href="http://projects.gnome.org/evince/">Evince</a>', ...
      'Matlab figure to be opened with Matlab or <a href="http://ifit.mccode.org">iFit</a>. Use <i>set(gcf,''visible'',''on'')</i> after loading.', ...
      'Portable Document File to be viewed with <a href="http://get.adobe.com/fr/reader/">Acrobat Reader</a> or <a href="http://projects.gnome.org/evince/">Evince</a>.' , ...
    '<a href="https://en.wikipedia.org/wiki/Encapsulated_PostScript">Encapsulated PostScript</a>, to be viewed with e.g. <a href="http://inkscape.org/">Inkscape</a>, <a href="http://projects.gnome.org/evince/">Evince</a>.', ...
    '<a href="https://fr.wikipedia.org/wiki/Scalable_Vector_Graphics">Scalable Vector Graphics</a> image, to be viewed with Chrome/Firefox, <a href="http://inkscape.org/">Inkscape</a>, <a href="http://www.gimp.org/>GIMP.</a>, <a href="http://projects.gnome.org/evince/">Evince</a>' ];
    
    % add image and links to exported Plot formats
    if ~isempty(dir([ basename '.png' ]))
      fprintf(fid, '<div style="text-align: center;"><a href="%s"><img src="%s" align="middle"></a><br>\n<i>Model: %s</i><br></div>\n', ...
        [ basename_img '.png' ], ...
        [ basename_img '.png' ], titl);
    end
    
    % display list of available Model formats, as well as suggested software to use
   fprintf(fid, '<h4>Model Exported to: </h4><ul>\n');
    for index=1:numel(export)
      if ~isempty(dir([ basename '.' export{index} ]))
        fprintf(fid, [ '<li><b><a href="' basename_img '.' export{index} '">' export{index} '</a></b>: ' ...
          export_label{index} '</li>\n' ]);
      end
    end
    fprintf(fid, '</ul></p>\n');
  
  case 'expression'
  
    % Model expression (details)
    fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', Section, Section);
    
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
    fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', Section, Section);
    
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
  end % switch
  
  if ~isempty(message), publish(self, filename, 'message', message); end
  fclose(fid);


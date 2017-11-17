function filename = publish(self, filename)
% iFunc_McCode: publish: generate a readable HTML document and export results to
% many file formats
%
%  publish(s)
%    generates a readable report for a 4D S(q,w) object
%  publish(s,'force')
%    same as above, but make sure any existing report is overwritten
%  publish(s,filename)
%    explicitly write the report into given file


% section defaults to {'header','crystal','calculator','properties',
%   'plot','thermochemistry','download','footer'};
% other possible sections are: 'status', 'error','optimize'

  
  section_all  = {'header', 'toc', 'parameters', ...
      'plot','model', ...
      'monitors', ...
      'download', ...
      'footer'}; 

  if nargin < 2, filename = []; end
  if nargin < 3, section  = []; end
  if nargin < 4, message  = []; end
  
  if strcmp(filename,'force') && isempty(section)
    filename=''; section='force';
  end
  
  UserData = self.UserData; % this is where we search for information
  if isfield(UserData, 'options'), options = UserData.options; else options = []; end
  if isempty(filename) && isfield(options, 'report') && ~isempty(options.report)
    filename = options.report;
  end
  % overwrite when in 'force' mode
  if any(strcmp(section, 'force'))
    if ~isempty(dir(filename)), delete(filename); end
    if iscell(section), section(strcmp(section,'force')) = [];
    else section = []; end
  end
  if isempty(section), section = section_all; end
  
  if ~iscell(section), section = { section }; end
  
  % write content ==============================================================
  
  for index=1:numel(section)
    % open the file
    try
      fid = fopen(filename, 'a+');  % create or append to file
    catch
      fid = -1;
    end
  
    switch strtok(lower(section{index}))
    case {'header','model','download','footer','message','force'}
      % these are supported by the iFunc/publish method
      %   this should e.g. be the first call with section='header'
      %   which creates the file and its 'img' directory
      filename = publish@iFunc(self, filename, section{index});
      
    case 'table'
      if ischar(message)
        dummy = publish@iFunc(self, filename, 'table', message);
        message = [];
      end
      
    case 'plot'
      % Model plot ***************************************************************
      fprintf(fid,'<h2><a name="%s"></a>%s</h2>\n', 'plot', 'Plot');
      
      titl   = self.Name;
      titl(titl=='<')='[';
      titl(titl=='>')=']';
      desc   = self.Description;
      desc(desc=='<')='[';
      desc(desc=='>')=']';
      basename     = fullfile(options.target, 'img', [ 'iFit_Model_' self.Tag ]);
      
      f = figure('Visible','off', 'Name', [ 'iFit_Model_' self.Tag ]);
      h = plot(self); axis tight;

      % create output from the figure: png pdf fig
      set(f,'Renderer','zbuffer');
      saveas(f, basename, 'fig');
      saveas(f, basename, 'png');
      saveas(f, basename, 'pdf');
      saveas(f, basename, 'epsc');
      plot2svg([ basename '.svg' ], f, 'jpg');
      
      figure2xhtml([ basename '.x3d' ], f, struct('interactive',true, ...
        'output', 'x3d','title',titl,'Description',desc));
      figure2xhtml([ basename '.xhtml' ], f, struct('interactive',true, ...
        'output', 'xhtml','title',titl,'Description',desc));
      close(f);
        
    case 'parameters'
      % call the default publish (parameters) but add the Parameters_Constant
      dummy = publish@iFunc(self, filename, 'parameters');
      
      if isfield(self.UserData,'Parameters_Constant') && ~isempty(self.UserData.Parameters_Constant)
        dummy = publish@iFunc(self, filename, 'message', '<p>Additional instrument parameters</p><br>');
        publish_properties(self, fid, self.UserData.Parameters_Constant);
      end
      dummy = publish@iFunc(self, filename, 'message', '<p>Simulation parameters</p><br>');
      publish_properties(self, fid, options);
      
      % save the instrument description into the target directory
      UD = self.UserData;
      if ~isempty(UD.instrument_source) && ...
        isempty(dir(fullfile(options.target, 'img', options.instrument)))
        f=fopen(fullfile(options.target, 'img', options.instrument), 'w');
        if f==-1, return; end
        fprintf(f, '%s\n', UD.instrument_source);
        fclose(f);
      end
      % publish it using the publish native call
      html = publish(fullfile(options.target, 'img', options.instrument), ...
        struct('evalCode',false,'outputDir',fullfile(options.target, 'img')));
      % add a link to the instrument source code
      [p,f,e] = fileparts(html);
      dummy = publish@iFunc(self, filename, 'message', ...
          sprintf([ '<br>Instrument source code <a href="%s">%s</a>' ...
            ' [<a href="%s">raw description file</a>].\n ' ], ...
          fullfile('img', [f e]), options.instrument, fullfile('img', options.instrument)));
      
    case 'toc'
      
      % add the list of tags
      publish_toc(section, fid, options);

    case 'status'
    
      % display the current status of the computation (start, ETA, percentage...)
      if isfield(options, 'status')
        if ~ischar(options.status) options.status = num2str(options.status); end
        dummy = publish@iFunc(self, filename, 'message', ...
          sprintf('<br>[%s] %s %s\n', datestr(now), options.status));
      end
    
    case 'monitors'
    
      self = publish_monitors(self, fid, options);
      dummy = publish@iFunc(self, '', 'table', 'McCode_OverView');
      
      % list the simulated data 'monitors'
      h5_files = saveas(self.UserData.monitors, fullfile(options.target, 'img', 'McCode_OverView.h5'), 'mantid');
      if ischar(h5_files), h5_files = { h5_files }; end
      mon_all = self.UserData.monitors;
      
      dummy = publish@iFunc(self, filename, 'message', '<p><b>Monitor files:</b></p><ul>' );
      for index=1:numel(h5_files)
        if numel(mon_all) == 1, this_mon = mon_all; else this_mon = mon_all(index); end
        [p,f,e]    = fileparts(h5_files{index});
        [p1,f1,e1] = fileparts(get(this_mon, 'Source'));
        dummy = publish@iFunc(self, filename, 'message', ...
          sprintf('<li><a href="%s">%s</a> from <b>%s</b> [HDF5 Mantid NeXus]</li>\n', ...
            fullfile('img', [f e ]), [f e], [f1 e1] ));
      end
      if ~isempty(dir(fullfile(options.target, 'sim')))
        dummy = publish@iFunc(self, filename, 'message', '</ul><p>Raw monitor files are available in: <a href="sim">sim</a></p>' );
      else
        dummy = publish@iFunc(self, filename, 'message', '</ul><br>' );
      end  
      
    case 'error'
      % display error message
      publish(self, filename, 'status', ...
        sprintf('<hr><h2>ERROR: %s %s FAILED</h2>\n', calculator, configuration) );
    otherwise
      disp([ mfilename ': INFO: ignoring section ' section{index} ]);
    end % switch
    if fid ~= -1, fclose(fid);
    if ~isempty(message), dummy = publish@iFunc(self, filename, 'message', message); end
    
    % only one custom message per call
    if numel(section) > 1 && ~isempty(message), message = []; end
    
    % update filename into current object
    if ~isempty(filename)
      options.report   = filename; end
      options.target   = fileparts(filename);
      UserData.options = options;
      self.UserData    = UserData;
    end
    
    if strcmp(strtok(lower(section{index})),'error')
      return; % abort on error
    end
    disp([ mfilename ': Generated section ' filename '#' strtok(lower(section{index})) ]);
    
  end % for (section)
  
  if ~isempty(inputname(1))
    assignin('caller',inputname(1),self); % update in original object
  end
  
  % last, open web browser when no output requested ------------------------------
  if nargout == 0
    webbrowser(filename, 'system');
  end

% ==============================================================================
function publish_properties(self, fid, options)

  if fid == -1, return; end
  
  % display some more information about the model
  toadd = fieldnames(options);
  fprintf(fid, '<table style="text-align: left; width: 80%%;" border="1" cellpadding="2" cellspacing="2">\n');
  for index=1:numel(toadd)
    this = options.(toadd{index});
    if isnumeric(this) && ndims(this) <= 2
      this = double(this);
      if numel(this) > 100
        this1 = this(1:70);
        this2 = this((end-20):end);
        fprintf(fid, '  <tr><td><b>%s</b></td><td>%s ... %s</td></tr>\n', ...
          toadd{index}, mat2str(this1), mat2str(this2));
      else
        fprintf(fid, '  <tr><td><b>%s</b></td><td>%s</td></tr>\n', toadd{index}, mat2str(this));
      end
    elseif ischar(this)
      fprintf(fid, '  <tr><td><b>%s</b></td><td>%s</td></tr>\n', toadd{index}, this');
    end
  end
  fprintf(fid, '</table></p>\n');

  
% ==============================================================================
function self = publish_monitors(self, fid, options)
  % the S(0kl,w)
  
  if fid == -1, return; end
  
  f = figure('Visible','off');
  h = subplot(self, 'view2 tight');
  % export overview plot
  % create output from the figure: png pdf fig
  basename = fullfile(options.target, 'img', 'McCode_OverView');
  set(f,'Renderer','zbuffer');
  saveas(f, basename, 'fig');
  saveas(f, basename, 'png');
  saveas(f, basename, 'pdf');
  saveas(f, basename, 'epsc');
  plot2svg([ basename '.svg' ], f, 'jpg');
  close(f);

  fprintf(fid, '<h3><a name="monitors"></a>The simulation results (monitors)</h3>\n');
  fprintf(fid, '<p>Here is a view of the generated monitors using the <a href="#parameters">above</a> parameters.<br>\n');
    
  McCode_OverView = self.UserData.monitors;
  builtin('save',    [ basename '.mat' ], 'McCode_OverView');

% ==============================================================================
function publish_toc(section, fid, options)
  % table of contents
  if fid == -1, return; end
  if ~iscellstr(section) || numel(section) <= 1, return; end
  
  fprintf(fid, '<hr>Table of contents:<br><ul>\n');
  for index=1:numel(section)
    if any(strcmp(section{index}, {'header','footer','toc','plot'})), continue; end
    fprintf(fid, '<li><a href="#%s">%s</a></li>\n', section{index}, section{index});
  end
  fprintf(fid, '</ul><hr>\n');



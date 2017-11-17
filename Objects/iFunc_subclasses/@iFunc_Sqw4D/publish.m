function filename = publish(self, filename, section, message)
% iFunc_Sqw4D: publish: generate a readable HTML document and export results to
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
      'crystal', ...
      'calculator', ...
      'optimize', ...
      'properties', ...
      'plot','model', ...
      'thermochemistry', ...
      'plot3', ...
      'powder', ...
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
    case {'header','plot','model','parameters','download','footer','message','force'}
      % these are supported by the iFunc/publish method
      %   this should e.g. be the first call with section='header'
      %   which creates the file and its 'img' directory
      filename = publish@iFunc(self, filename, section{index});
    case 'table'
      if ischar(message)
        dummy = publish@iFunc(self, filename, 'table', message);
        message = [];
      end
    case 'toc'
      
      % add the list of tags
      publish_toc(section, fid, options);
      
    case 'crystal'
    
      publish_crystal(self, fid, options);
    
    case 'calculator'
    
      publish_calculator(self, fid, options);
      
    case 'optimize'
    
      publish_optimize(self, fid, options);
    
    case 'status'
    
      % display the current status of the computation (start, ETA, percentage...)
      if isfield(options, 'status')
        if ~ischar(options.status) options.status = num2str(options.status); end
        dummy = publish@iFunc(self, filename, 'message', ...
          sprintf('<br>[%s] %s %s\n', datestr(now), options.status));
      end
    
    case 'properties'
    
      publish_properties(self, fid, options);
    
    case 'dos','thermochemistry'
    
      publish_dos(self, fid, options);
    
    case 'plot3'
    
      self = publish_eval_3D(self, fid, options);
    
    case 'powder'
    
      publish_eval_powder(self, fid, options)
      
    case 'error'
      % display error message
      if isfield(options, 'configuration') configuration = options.configuration;
      else configuration=self.Name; end
      if isfield(options, 'calculator')    calculator = options.calculator;
      else calculator=[]; end
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
  if nargout == 0 && numel(section) > 1
    webbrowser(filename, 'system');
  end

  
% ==============================================================================
function publish_crystal(self, fid, options)
  % introduction, general information
  
  if fid == -1, return; end

  % introduction: All S(q,w) 4D objects
  message = { [ '<p>This is an ' class(self) ' object describing a 4D S(q,w) <a href="https://en.wikipedia.org/wiki/Dynamic_structure_factor">dynamic structure factor</a>. ' ], ...
    'It represents the double Fourier transform of the space-time pair correlation distribution function into the reciprocal space-energy frame. ', ...
    'In practice the S(q,w=0) is the structure factor, where as the w<0 and w>0 parts are the inelastic scattering contributions.</p>' };
  fprintf(fid, '%s<br>\n', message{:});
    
  % test if this is an ASE/PhonoPy derived object
  if isfield(options, 'use_phonopy')
    % sqw_phonons objects
    fprintf(fid, '<p>This page presents the results of the estimate of phonon dispersions (lattice dynamics) in a single crystal, using a DFT code (the "calculator"). From the initial atomic configuration (geometry), each atom in the lattice cell is displaced by a small quantity. The displaced atom then sustains a, so called Hellmann-Feynman, restoring force to come back to the stable structure. The dynamical matrix is obtained from these forces, and its eigen-values are the energies of the vibrational modes in the crystal.</p>\n');
    fprintf(fid, '<p>This computational resource is provided by <a href="http://ifit.mccode.org">iFit</a>, with the <a href="http://ifit.mccode.org/Models_Phonons.html"><b>Phonons</b> Model</a>, which itself makes use of the <a href="https://wiki.fysik.dtu.dk/ase">Atomic Simulation Environment (ASE)</a>.\n');
    if options.use_phonopy
      fprintf(fid, ' In addition, the <a href="https://atztogo.github.io/phonopy/">PhonoPy</a> package is used to compute force constants (faster, more accurate).\n');
    end
    fprintf(fid, '<p>This report summarizes the initial crystal geometry and the calculator configuration. When the simulation ends successfully, the lower part presents the S(hkl,w) dispersion curves as plots and data sets, the model used (as a Matlab object), and the density of states. These results correspond to the coherent inelastic part of the dynamic structure factor S(hkl,w), for vibrational modes in the harmonic approximation. In the following, we use energy unit in meV = 241.8 GHz = 11.604 K = 0.0965 kJ/mol, and momentum is given in reduced lattice units.</p>\n');
    fprintf(fid, '<p><b>Limitations:</b> The accuracy of the model depends on the parameters used for the computation, e.g. energy cut-off, k-points grid, smearing, ...</p>\n');
  end
  
  % append crystal information
  if isfield(options, 'configuration')
    fprintf(fid, '<h2><a name="crystal"></a>Atom/molecule initial structure</h2>\n');
    if ~isempty(dir(options.configuration))
      [p,f,e] = fileparts(options.configuration);
      if isfield(options, 'target') && isempty(dir(fullfile(options.target, [f e])))
        copyfile(options.configuration, fullfile(options.target, [f e])); 
      end
      fprintf(fid, '<p>In this section, we present the crystallographic lattice cell used for the computation, both as plots and structure files for use with other software.</p>\n');
      fprintf(fid, 'Initial description: <a href="%s">%s</a><br>\n', [f e], [f e]);
    elseif ischar(options.configuration)
      fprintf(fid, 'Initial description: %s\n', options.configuration);
    end
    fprintf(fid, '<br>\n');
    % look for exported 'configuration' (structure) files
    if isfield(options, 'target')
      publish(self, '', 'table', 'configuration');
      
      if ~isempty(dir(fullfile(options.target, 'sqw_phonons_check.py')))
        fprintf(fid, 'The python/ASE script used to import the initial material description is:\n<br>');
        fprintf(fid, '<ul><li><a href="sqw_phonons_check.py">sqw_phonons_check.py</a></li></ul>\n');
      end
    end
  end
  
% ==============================================================================
function publish_calculator(self, fid, options)
  
  % tag=calc: append calculator configuration
  
  if fid == -1, return; end

  if isfield(options, 'calculator')
    fprintf(fid, '<h2><a name="calculator"></a>Calculator configuration</h2>\n');
    % general introduction
    fprintf(fid, 'This is a quick desription of the calculator used to compute the lattice energy and forces.<br>\n');
    if isfield(options.link) && isstruct(options.link)
      fprintf(fid, 'We are using the calculator: <a href="%s"><b>%s</b></a><br>\n', ...
        options.link.(options.calculator), upper(options.calculator));
    else
      fprintf(fid, 'We are using the calculator: <b>%s</b><br>\n', ...
        upper(options.calculator));
    end
    
    % clean 'options' from empty and 0 members, scripts, etc
    op           = options;
    toremove = {'available','cite','autoplot','available','gui','htmlreport','link'};
    for index=1:numel(toremove)
      if isfield(op, toremove{index})
        op           = rmfield(op, toremove{index});
      end
    end
    for f=fieldnames(op)'
      if isempty(op.(f{1})) ...
        || (isnumeric(op.(f{1})) && isscalar(op.(f{1})) && (op.(f{1}) == 0 || ~isfinite(op.(f{1})))) ...
        || strncmp(f{1}, 'script',6)
        op = rmfield(op, f{1});
      end
    end
  
    if ~isempty(op)
      fprintf(fid, [ 'The calculator configuration:<br>' ...
        '<div style="margin-left: 40px;"><pre>%s</pre></div>\n' ], ...
        class2str(' ',op));
    end
    fprintf(fid, '<hr>\n');
  end
  
% ==============================================================================
function publish_optimize(self, fid, options)

  if fid == -1, return; end
  % display result of the optimization
  if ~isfield(options, 'optimizer') || isempty(options.optimizer), return; end
  if ~isfield(options, 'target') || isempty(dir(options.target)), return; end
  
  fprintf(fid, '<h2><a name="optimize"></a>Optimized atom/molecule configuration</h2>\n');
  fprintf(fid, '<p>In this section, we present the results of the optimization procedure, of the initial structure. The lattice parameters are kept constant. The optimized structure is used further for the force estimate.</p>\n');
  fprintf(fid, '<br>\n');
  publish(self, '', 'table', 'optimized');
  
  if ~isempty(dir(fullfile(options.target, 'sqw_phonons_optimize.py')))
    fprintf(fid, 'The python/ASE script used to optimize the initial material structure is\n');
    fprintf(fid, '<a href="sqw_phonons_optimize.py">sqw_phonons_optimize.py</a>\n');
  end
  fprintf(fid, '<hr>\n');
  
  
% ==============================================================================
function publish_properties(self, fid, options)

  if fid == -1, return; end
  
  % display some information about the system (energy, structure, etc...)
  if isfield(self.UserData, 'properties')
    fprintf(fid, '<h2><a name="properties">Physical properties</h2>\n');
    properties = self.UserData.properties;
    toadd = fieldnames(properties);
    fprintf(fid, '<table style="text-align: left; width: 80%%;" border="1" cellpadding="2" cellspacing="2">\n');
    for index=1:numel(toadd)
      this = properties.(toadd{index});
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
  end
  
  if ~isempty(dir(fullfile(options.target, 'sqw_phonons_forces_iterate.py')))
    fprintf(fid, 'The python/ASE scripts used to compute the forces are:\n<br>');
    fprintf(fid, '<ul><li><a href="sqw_phonons_eval.py">sqw_phonons_forces_iterate.py</a></li>\n');
    fprintf(fid, '<li><a href="sqw_phonons_eval.py">sqw_phonons_forces_finalize.py</a></li></ul>\n');
  end
  
% ==============================================================================
function publish_dos(self, fid, options)
  % display the vDOS. Model must have been evaluated once to compute DOS
  
  if fid == -1, return; end
  try
    Sqw4D_DOS = dos(self);
  catch ME
    disp([ mfilename ': ERROR: can not get the vDOS for object ' class(self) ]);
    return
  end
  try
    [thermo, fig] = thermochemistry(self, 1:1000, 'plot');
  catch ME
    disp([ mfilename ': WARNING: can not get the Thermo-Chemistry quantities for object ' class(self) ]);
    thermo = [];
  end
  if ~isempty(Sqw4D_DOS)
    if isfield(thermo, 'entropy'),          Sqw4D_Thermo_S = thermo.entropy; 
    else                                    Sqw4D_Thermo_S=[]; end
    if isfield(thermo, 'internal_energy'),  Sqw4D_Thermo_U = thermo.internal_energy;  
  else                                      Sqw4D_Thermo_U=[]; end
    if isfield(thermo, 'helmholtz_energy'), Sqw4D_Thermo_F = thermo.helmholtz_energy;  
    else                                    Sqw4D_Thermo_F=[]; end
    if isfield(thermo, 'heat_capacity'),    Sqw4D_Thermo_Cv= thermo.heat_capacity;  
    else                                    Sqw4D_Thermo_Cv=[]; end
    
    fprintf(fid, '<h3><a name="thermochemistry"></a>The vibrational density of states (vDOS)</h3>\n');
    fprintf(fid, '<p>The <a href="https://en.wikipedia.org/wiki/Density_of_states">vibrational density of states</a> (aka phonon spectrum) can be defined as the velocity auto-correlation function (VACF) of the particles. ');
    if ~isempty(Sqw4D_Thermo_S)
      fprintf(fid, 'The <a href="https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html#background">Thermodynamic quantities</a> have also been computed.');
    end
    fprintf(fid, '</p>\n');
    if ishandle(fig)
      set(fig,'Renderer','zbuffer')
      saveas(fig, fullfile(options.target,  'img', 'Sqw4D_ThermoChem.fig'),  'fig');
      set(fig,'Visible','off');
      plot2svg(fullfile(options.target,     'img', 'Sqw4D_ThermoChem.svg'), fig, 'jpg');
      print(fig,  fullfile(options.target,  'img', 'Sqw4D_ThermoChem.png'),'-dpng');
      print(fig,  fullfile(options.target,  'img', 'Sqw4D_ThermoChem.pdf'),'-dpdf');
      
      publish(self, '', 'table', 'Sqw4D_ThermoChem');
      close(fig); fig=[];
    end
    
    properties = {'Sqw4D_DOS','Sqw4D_Thermo_S','Sqw4D_Thermo_U','Sqw4D_Thermo_F','Sqw4D_Thermo_Cv'};
    for index=1:numel(properties)
      name = properties{index};
      if isempty(name), continue; end
      try
        val = eval(name); % these should be iData objects
        builtin('save', fullfile(options.target,  'img', [ name '.mat']), name);
        save(val, fullfile(options.target, 'img', [ name '.dat']), 'dat data');
        save(val, fullfile(options.target, 'img', [ name '.h5']), 'mantid');
      catch
        disp([ mfilename ': ERROR when exporting ' name ' into HTML. Skipping.' ])
        continue
      end
      switch name
      case 'Sqw4D_DOS'
        fprintf(fid, 'The vibrational spectrum is available below:</p><br>\n');
      case 'Sqw4D_Thermo_S'
        fprintf(fid, 'The entropy S=-dF/dT [J/K/mol] is available below:</p><br>\n');
      case 'Sqw4D_Thermo_U'
        fprintf(fid, 'The internal energy U [J/mol] is available below:</p><br>\n');
      case 'Sqw4D_Thermo_F'
        fprintf(fid, 'The Helmholtz free energy F=U-TS=-kT lnZ [J/mol] is available below:</p><br>\n');
      case 'Sqw4D_Thermo_Cv'
        fprintf(fid, 'The molar specific heat at constant volume Cv=dU/dT [J/K/mol] is available below:</p><br>\n');
      end
      publish(self, '', 'table', name);
    
    end
    fprintf(fid, '<hr>\n');
  end % if Sqw4D_DOS
  if ishandle(fig), close(fig); end
  
% ==============================================================================
function self = publish_eval_3D(self, fid, options)
  % the S(0kl,w)
  
  if fid == -1, return; end
  
  maxFreq = max(self);
  % evaluate the 4D model onto a 3D mesh filling the Brillouin zone [0:0.5 ] at QH=0
  qk=linspace(0,1.5,50); qh=0; ql=qk; 
  w =linspace(0.01,maxFreq*1.2,51);
  Sqw4D_0KLE =iData(self,[],qh,qk,ql',w); Sqw4D_0KLE = squeeze(Sqw4D_0KLE);

  fprintf(fid, '<h3><a name="plot3"></a>The dispersion in 3D S(h=0,k,l,w)</h3>\n');
  fprintf(fid, '<p>In order to view this 4D data set, we represent it on the QH~0 plane as a 3D volume data set. The intensity level is set as log10[S(QH=%g,QK,QL,w)].<br>\n', qh(1));
    
  fprintf(fid, '<ul><li>qh=%g (QH in rlu)</li>\n', qh(1)); % axis1
  fprintf(fid, '<li>qk=[%g:%g] with %i values (QK in rlu)</li>\n', xlim(Sqw4D_0KLE), size(Sqw4D_0KLE,2));     % axis2
  fprintf(fid, '<li>ql=[%g:%g] with %i values (QL in rlu)</li>\n', ylim(Sqw4D_0KLE), size(Sqw4D_0KLE,3));
  fprintf(fid, '<li>w=[%g:%g] with %i values (energy in meV, vertical)</li></ul></p>\n', zlim(Sqw4D_0KLE), size(Sqw4D_0KLE,4));
    
    fprintf(fid, '<p>The QH=%g data set is available in the folowing formats (log10 of the data set except for DAT, HDF5 and MAT files):<br>\n', qh(1));
    
  builtin('save',    fullfile(options.target, 'img', 'Sqw4D_0KLE.mat'), 'Sqw4D_0KLE');
  saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.dat'), 'dat data');
  saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.h5'), 'mantid');
  
  Sqw4D_0KLE = log10(Sqw4D_0KLE);
  saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.png'), 'png', 'tight');
  saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.fig'), 'fig', 'plot3 tight');
  saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.vtk'));
  saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.mrc'));
  
  % the PDF export may crash in deployed version
  if ~isdeployed
    saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.pdf'), 'pdf', 'tight');
  end
  
  % modify aspect ratio to fit in a cube for X3D/XHTML
  saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.xhtml'), 'xhtml','axes auto');
  saveas(Sqw4D_0KLE, fullfile(options.target, 'img', 'Sqw4D_0KLE.x3d'), 'x3d','axes auto');
  
  publish(self, '', 'table', 'Sqw4D_0KLE');
  
  if ~isempty(dir(fullfile(options.target, 'sqw_phonons_eval.py')))
    fprintf(fid, 'The python/ASE script used to evaluate the HKLE dispersion is:\n<br>');
    fprintf(fid, '<ul><li><a href="sqw_phonons_eval.py">sqw_phonons_eval.py</a></li></ul>\n');
  end
 
% ==============================================================================
function Phonon_powder = publish_eval_powder(self, fid, options)

  if fid == -1, return; end
  % the S(hklw) radial average
  Sqw2D_Powder = powder(self); % a 2D iFunc
  Sqw2D_Powder_eval = iData(Sqw2D_Powder,[],linspace(0,4,30),linspace(0,max(self)*1.2,51));
  log_Sqw2D_Powder = log10(Sqw2D_Powder_eval);
  
  fprintf(fid, '<h3><a name="powder"></a>The powder average S(q,w)</h3>\n');
  fprintf(fid, 'The powder average S(q,w) is shown below:<br>\n');
    
  builtin('save', fullfile(options.target, 'img', 'Sqw2D_Powder.mat'), 'Sqw2D_Powder');
  saveas(Sqw2D_Powder_eval, fullfile(options.target, 'img', 'Sqw2D_Powder.dat'), 'dat data');
  saveas(Sqw2D_Powder_eval, fullfile(options.target, 'img', 'Sqw2D_Powder.h5'), 'mantid');
  
  saveas(log_Sqw2D_Powder, fullfile(options.target, 'img', 'Sqw2D_Powder.png'),'png data');
  saveas(log_Sqw2D_Powder, fullfile(options.target, 'img', 'Sqw2D_Powder.fits'),'fits');
  saveas(log_Sqw2D_Powder, fullfile(options.target, 'img', 'Sqw2D_Powder.tiff'),'tiff');
  saveas(log_Sqw2D_Powder, fullfile(options.target, 'img', 'Sqw2D_Powder.fig'), 'fig', 'tight');
  saveas(log_Sqw2D_Powder, fullfile(options.target, 'img', 'Sqw2D_Powder.pdf'), 'pdf', 'tight view2');
  saveas(log_Sqw2D_Powder, fullfile(options.target, 'img', 'Sqw2D_Powder.svg'), 'svg', 'view2 tight');
  saveas(log_Sqw2D_Powder, fullfile(options.target, 'img', 'Sqw2D_Powder.x3d'), 'x3d', 'axes auto');
  
  publish(self, '', 'table','Sqw2D_Powder');
  

% ==============================================================================
function publish_toc(section, fid, options)
  % table of contents
  if fid == -1, return; end
  if ~iscellstr(section) || numel(section) <= 1, return; end
  
  fprintf(fid, '<hr>Table of contents:<br><ul>\n');
  for index=1:numel(section)
    if any(strcmp(section{index}, {'header','footer','toc','plot'})), continue; end
    if ~strcmp(section{index},'optimize') || ...
       ~strcmp(section{index},'optimize') || ...
       (strcmp(section{index},'optimize') && isfield(options, 'optimizer') && ...
       ~isempty(options.optimizer)) || ...
       (strcmp(section{index},'calculator') && isfield(options, 'calculator') && ...
       ~isempty(options.calculator))
      fprintf(fid, '<li><a href="#%s">%s</a></li>\n', section{index}, section{index});
    end
  end
  fprintf(fid, '</ul><hr>\n');


% ==============================================================================


  


function filename = sqw_phonons_htmlreport(filename, step, options, data)
% write an HTML report on a sqw_phonons computation
%
% sections:
% * calculator and system informations (from configuration, exported files, options)
% * optimisation status
% * on-going calculation status (when possible, e.g. QE/PHON will write ETA)
% * status at end of forces: success/failure, properties
% * success:
%     estimate max spectrum energy with coarse grid.        write Emax to report
%     compute kpath using the spacegroup -> path            export kpath iData
%     evaluate model on finer grid with [0 Emax]            export 4D dat
%                                                           export 3D H=0 iData
%     compute powder spectrum                               export powder iData
%     compute vDOS                                          export vDOS iData
%     model                                                 export YAML/MAT/JSON
if isempty(options.htmlreport) || ~options.htmlreport
  return
end
if isempty(step) || ~ischar(step), return; end

if isempty(filename)
  filename = fullfile(options.target, 'sqw_phonons.html');
end

% open the filename
if strcmp(step, 'init')
  fid = fopen(filename, 'w');
else
  fid = fopen(filename, 'a');
end
if fid == -1, return; end
options.report = filename;

switch lower(step)
case 'create_atoms'
  % display result of initial material import/export
  sqw_phonons_htmlreport_create_atoms(fid, options);
case 'optimize'
  % display result of the optimization
  sqw_phonons_htmlreport_optimize(fid, options);
case 'status'
  % display status of the computation (starting, ETA, percentage, ...)
  sqw_phonons_htmlreport_status(fid, options);
case 'results'
  %
  if isdeployed || ~usejava('jvm') || ~usejava('desktop')
    disp([ 'sqw_phonons: generating HTML report in ' options.report ]);
  else
    disp([ 'sqw_phonons: generating HTML report in <a href="' options.report '">' options.report '</a>' ]);
  end
  Phonon_Model = sqw_phonons_htmlreport_model(fid, options);
  if isempty(Phonon_Model), return; end
  [maxFreq, Phonon_Model]= sqw_phonons_htmlreport_max_spectrum(fid, options, Phonon_Model);
  [grid4D,  Phonon_Model]= sqw_phonons_htmlreport_eval_4D(fid, options, Phonon_Model, maxFreq);
  try
                 sqw_phonons_htmlreport_kpath(fid, options, Phonon_Model, maxFreq);
  catch ME
                 disp(getReport(ME))
  end
                 sqw_phonons_htmlreport_dos(fid, options, Phonon_Model);
                 sqw_phonons_htmlreport_eval_3D(fid, options, grid4D);
  sqw_phonons_htmlreport_eval_powder(fid, options, grid4D, maxFreq);
  
case 'download'
  sqw_phonons_htmlreport_download(fid, options)
case 'error'
  sqw_phonons_htmlreport_error(fid, options, data)
otherwise
  % pass
end
fclose(fid);

% ==============================================================================
function sqw_phonons_htmlreport_create_atoms(fid, options)
  
  % get icons, and clean up 'options' for display
  [logo, link, op] = sqw_phonons_htmlreport_init(options);
  
  % open the report HTML page and write header, title, date, ...
  fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n');
  fprintf(fid, '<html>\n<head>\n<title>%s: %s with %s</title>\n', ...
    'sqw_phonons', options.configuration, options.calculator);
  fprintf(fid, '</head>\n');
  fprintf(fid, '<body><div style="text-align: center;">\n');
  fprintf(fid, '<a href="http://ifit.mccode.org"><img title="ifit.mccode.org" src="iFit-logo.png" align="middle" height=100></a>\n');
  fprintf(fid, '<a href="https://wiki.fysik.dtu.dk/ase" title="wiki.fysik.dtu.dk/ase"><img src="ase256.png" align="middle" height=100></a></br>\n');
  fprintf(fid, '<h1>%s: %s with %s</h1></div>\n', 'sqw_phonons', options.configuration, options.calculator);
  % general introduction
  fprintf(fid, 'Start Date: %s.<br>\nComputer: %s.<br>\n', datestr(now), computer);
  fprintf(fid, 'Stored in: <a href="%s">%s</a><br>\n', options.target, options.target);
  
  % introduction
  fprintf(fid, '<p>This page presents the results of the estimate of phonon dispersions (lattice dynamics) in a single crystal, using a DFT code (the "calculator"). From the initial atomic configuration (geometry), each atom in the lattice cell is displaced by a small quantity. The displaced atom then sustains a, so called Hellmann-Feynman, restoring force to come back to the stable structure. The dynamical matrix is obtained from these forces, and its eigen-values are the energies of the vibrational modes in the crystal.</p>\n');
  fprintf(fid, '<p>This computational resource is provided by <a href="http://ifit.mccode.org">iFit</a>, with the <a href="http://ifit.mccode.org/Models.html#mozTocId990577"<b>sqw_phonon</b> Model</a>, which itself makes use of the <a href="https://wiki.fysik.dtu.dk/ase">Atomic Simulation Environment (ASE)</a>.\n');
  if options.use_phonopy
    fprintf(fid, ' In addition, the <a href="https://atztogo.github.io/phonopy/">PhonoPy</a> package is used to compute force constants (faster, more accurate).\n');
  end
  fprintf(fid, '<p>This report summarizes the initial crystal geometry and the calculator configuration. When the simulation ends successfully, the lower part presents the S(hkl,w) dispersion curves as plots and data sets, the model used (as a Matlab object), and the density of states. These results correspond to the coherent inelastic part of the dynamic structure factor S(hkl,w), for vibrational modes in the harmonic approximation. In the following, we use energy unit in meV = 241.8 GHz = 11.604 K = 0.0965 kJ/mol, and momentum is given in reduced lattice units.</p>\n');
  fprintf(fid, '<p><b>Limitations:</b> The accuracy of the model depends on the parameters used for the computation, e.g. energy cut-off, k-points grid, smearing, ...</p>\n');
  
  sqw_phonons_htmlreport_toc(fid, options);

  % append system configuration
  fprintf(fid, '<h2><a name="atom"></a>Atom/molecule configuration</h2>\n');
  fprintf(fid, '<p>In this section, we present the crystallographic lattice cell used for the computation, both as plots and structure files for use with other software.</p>\n');
  [p,f,e] = fileparts(options.configuration);
  fprintf(fid, 'Initial description: <a href="%s">%s</a><br>\n', [f e], [f e]);
  fprintf(fid, '<br>\n');
  sqw_phonons_htmlreport_table(fid, options, 'configuration');
  
  % tag=calc: append calculator configuration
  fprintf(fid, '<h2><a name="calc"></a>Calculator configuration</h2>\n');
  % general introduction
  fprintf(fid, 'This is a quick desription of the calculator used to compute the lattice energy and forces.<br>\n');
  if ~isempty(logo)
    fprintf(fid, 'We are using the calculator: <a href="%s"><img src="%s" height="80" align="middle"><b>%s</b></a> (<a href="%s">%s</a>).<br>\n', ...
      link, logo, upper(options.calculator), link, link);
  else
    fprintf(fid, 'We are using the calculator:  <b>%s</b><br>\n', upper(options.calculator));
  end
  
  fprintf(fid, 'The calculator configuration:<br><div style="margin-left: 40px;"><pre>%s</pre></div>\n', class2str(' ',op));
  if ~isempty(dir(fullfile(options.target, 'sqw_phonons_check.py')))
    fprintf(fid, 'The python/ASE script used to import the initial material description is:\n<br>');
    fprintf(fid, '<ul><li><a href="sqw_phonons_check.py">sqw_phonons_check.py</a></li></ul>\n');
  end
  fprintf(fid, '<hr>\n');
  
% ==============================================================================
function sqw_phonons_htmlreport_toc(fid, options)
  % table of contents
  fprintf(fid, '<hr>Table of contents:<br><ul>\n');
  fprintf(fid, '<li><a href="#atom">Crystallographic information</a></li>\n');
  fprintf(fid, '<li><a href="#calc">Calculator configuration</a></li>\n');
  if ~isempty(options.optimizer)
  fprintf(fid, '<li><a href="#optimize">The optimized structure</a></li>\n');
  end
  fprintf(fid, '<li><a href="#results">Results</a><ul>\n');
  fprintf(fid, '  <li><a href="#model">The Phonon Model S(hkl,w)</a></li>\n');
  fprintf(fid, '  <li><a href="#dos">The Phonon spectrum (vDOS)</a></li>\n');
  fprintf(fid, '  <li><a href="#kpath">The dispersion along principal directions</a></li>\n');
  fprintf(fid, '  <li><a href="#grid4d">The Model evaluated onto a 4D grid</a></li>\n');
  fprintf(fid, '  <li><a href="#grid3d">The Model evaluated onto a 3D grid (QH~0)</a></li>\n');
  fprintf(fid, '  <li><a href="#powder">The Powder average S(q,w)</a></li></ul></li>\n');
  fprintf(fid, '<li><a href="#zip">Download</a></li>\n');
  fprintf(fid, '</ul><hr>\n');

% ==============================================================================
function sqw_phonons_htmlreport_table(fid, options, name)
% write a Table which displays available exported 'configuration' files
  
  if isempty(name), name = 'configuration'; end
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
    if ~isempty(dir(fullfile(options.target,index1)))
      % the file exists
      if strcmp(index1, [ name '.png' ]) && ~isempty(dir(fullfile(options.target,[ name '.tiff' ])))
        fprintf(fid, '<div style="text-align:center"><a href="%s"><img src="%s" align="middle" title="%s"></a><br>(try the <a href="%s">TIFF file</a> in case the axes are not shown)</div><br>\n', index1, index1, index1, [ name '.tiff' ]);
      elseif strcmp(index1, [ name '.png' ])  % only the PNG, but not the TIFF
        fprintf(fid, '<div style="text-align:center"><a href="%s"><img src="%s" align="middle" title="%s"></a></div><br>\n', index1, index1, index1);
      elseif strcmp(index1, [ name '.html' ])
        % embed a frame
        fprintf(fid, '<div style="text-align:center"><iframe src="%s" align="middle" width="480" height="480"></iframe><br>%s<br>(<a href="%s" target=_blank>open in external window</a>)<br></div><br>\n', index1, index2, index1);
      elseif strcmp(index1, [ name '.xhtml' ])
        % embed an x3dom frame
        fprintf(fid, [ '<div style="text-align:center"><iframe src="%s" align="middle" width="700" height="850"></iframe><br>\n' ...
        '%s<br>(<a href="%s" target=_blank>open in external window</a>)<br>\n' ...
        'The <font color="blue">blue</font> axis is the Energy (meV), the <font color="red">red</font> axis is QK (rlu), the <font color="green">green</font> axis is QL (rlu).<br>\n' ...
        'You can rotate the model (left mouse button), zoom (right mouse button), and pan (middle mouse button).<br>\n' ], index1, index2, index1);
      else
        if flag_table == 0
          flag_table = 1;
          fprintf(fid, '<table style="text-align: left; width: 100%%;" border="1" cellpadding="2" cellspacing="2">\n');
        end
        fprintf(fid, '  <tr><td><a href="%s">%s</a></td><td>%s</td></tr>\n', index1, index1, index2);
      end
    end
  end
  if flag_table
    fprintf(fid, '</table><br>\n');
  end
  
% ==============================================================================
function [logo, link, op] = sqw_phonons_htmlreport_init(options)
  % retrieve a few external resources and clean up options
  
  copyfile(fullfile(ifitpath,'Docs','images','iFit-logo.png'), options.target);
  % http://ifit.mccode.org
  copyfile(fullfile(ifitpath,'Docs','images','ase256.png'), options.target);
  % https://wiki.fysik.dtu.dk/ase
  switch upper(options.calculator)
  case 'GPAW'
    logo='logo-gpaw.png';
    link='http://wiki.fysik.dtu.dk/gpaw';
    copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  case 'NWCHEM'
    logo='nwchem.png';
    link='http://www.nwchem-sw.org/';
    copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  case 'ELK'
    logo='elk.png';
    link='http://elk.sourceforge.net';
    copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  case {'DACAPO','JACAPO'}
    logo='jacapo.png';
    link='http://wiki.fysik.dtu.dk/dacapo';
    copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  case 'ABINIT'
    logo='abinit.png';
    link='http://www.abinit.org/';
    copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  case 'EMT'
    logo='emt.png';
    link='https://wiki.fysik.dtu.dk/ase/ase/calculators/emt.html';
    copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  case {'QUANTUMESPRESSO','QUANTUMESPRESSO_ASE'}
    logo='logo_qe.jpg';
    link='http://www.quantum-espresso.org/';
    copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  case 'VASP'
    logo='vasp.png';
    link='https://www.vasp.at/';
    copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  otherwise
    logo=''; link='';
  end
    
  % clean 'options' from empty and 0 members
  op           = options;
  op           = rmfield(op, 'available');
  op.htmlreport= 0;
  for f=fieldnames(op)'
    if isempty(op.(f{1})) ...
      || (isnumeric(op.(f{1})) && isscalar(op.(f{1})) && (op.(f{1}) == 0 || ~isfinite(op.(f{1})))) ...
      || strncmp(f{1}, 'script',6)
      op = rmfield(op, f{1});
    end
  end
  
% ==============================================================================
function sqw_phonons_htmlreport_optimize(fid, options)
  % display result of the optimization
  fprintf(fid, '<h2><a name="optimized"></a>Optimized atom/molecule configuration</h2>\n');
  fprintf(fid, '<p>In this section, we present the results of the optimization procedure, of the initial structure. The lattice parameters are kept constant. The optimized structure is used further for the force estimate.</p>\n');
  fprintf(fid, '<br>\n');
  sqw_phonons_htmlreport_table(fid, options, 'optimized');
  
  if ~isempty(dir(fullfile(options.target, 'sqw_phonons_optimize.py')))
    fprintf(fid, 'The python/ASE script used to optimize the initial material structure is\n');
    fprintf(fid, '<a href="sqw_phonons_optimize.py">sqw_phonons_optimize.py</a>\n');
  end
  fprintf(fid, '<hr>\n');
  
% ==============================================================================
function sqw_phonons_htmlreport_status(fid, options)
  % display the current status of the computation (start, ETA, percentage...)
  if isfield(options, 'status')
    if ~ischar(options.status), options.status = num2str(options.status); end
    fprintf(fid, '<br>[%s] %s\n', datestr(now), options.status);
  end

% ==============================================================================
function sqw_phonons_htmlreport_error(fid, options, message)
  % display error message
  fprintf(fid, '<hr><h2>ERROR: %s FAILED</h2>\n', options.configuration);
  fprintf(fid, '<p>%s</p>\n', message);

% ==============================================================================
function Phonon_Model = sqw_phonons_htmlreport_model(fid, options)
  % check if the Model has been created and export it.
  name = 'Phonon_Model';
  Phonon_Model = [];

  if isempty(dir(fullfile(options.target, [ name '.mat' ]))), return; end
  Phonon_Model = load(fullfile(options.target, [ name '.mat' ]));
  if isstruct(Phonon_Model)
    Phonon_Model = Phonon_Model.Phonon_Model;
  end
  
  if isempty(Phonon_Model) || ~isa(Phonon_Model, 'iFunc')
    sqw_phonons_htmlreport_error(fid, options, 'Model has not been generated.')
    return; 
  end
  
  % now we export the Model into other formats
  % indicate evaluated model, and print grid used
  fprintf(fid, '<h2><a name="results"></a>Computation completed</h2>\n');
  fprintf(fid, '<div style="text-align: center; color:#0000FF"><b>WELL DONE</b></div><br>');
  fprintf(fid, '<p>The computation has performed correctly. The forces and dynamical matrix of the lattice vibrations have been determined.</p>\n');
  
  fprintf(fid, '<p>The results are stored into a 4D <a href="http://ifit.mccode.org/iFunc.html">iFunc</a> object containing the dynamical matrix. This is a Matlab workspace (MAT-file).');
  fprintf(fid, '<ul><li><a href="%s">%s</a></li></ul>\n', ...
    'Phonon_Model.mat', 'Phonon_Model.mat');
  fprintf(fid, 'Load the Model under <a href="http://ifit.mccode.org">Matlab/iFit</a> with (this also works with the <a href="http://ifit.mccode.org/Install.html">standalone version of iFit</a> which does <b>not</b> require any Matlab license and installs on most systems): <ul><li>load(''<a href="%s">%s</a>'') <i>%% creates a "Phonon_Model" iFunc Model</i></li></ul>\n', 'Phonon_Model.mat', 'Phonon_Model.mat');
  fprintf(fid, 'Define axes for the evaluation grid in 4D, for instance:<ul><li>qh=linspace(0.01,.5,50); qk=qh; ql=qh; w=linspace(0.01,50,51);</li></ul>\n');
  fprintf(fid, 'Evaluate the Phonons as an <a href="http://ifit.mccode.org/iData.html">iData</a> object under Matlab/iFit with: <ul><li>s=iData(Phonons_Model, [], qh, qk, ql'', w) <i>%% evaluates the "Phonons" onto the grid, with default parameters, and return an iData object</i></li>\n');
  fprintf(fid,'<li>plot3(log(s(1,:,:,:))) <i>%% <a href="http://ifit.mccode.org/Plot.html">plot</a> the data set for QH=0.01 rlu in log scale</i></li>\n');
  fprintf(fid, '</ul></p>\n');
    
  fprintf(fid, '<p><a name="model"></a>The Phonon model object is available in the following formats.</p>\n');
  % generate the Model into additional formats
  try; save(Phonon_Model, fullfile(options.target,name), 'xml');  end
  try; save(Phonon_Model, fullfile(options.target,name), 'json'); end 
  try; save(Phonon_Model, fullfile(options.target,name), 'yaml'); end
  try; save(Phonon_Model, fullfile(options.target,name), 'm');    end
  sqw_phonons_htmlreport_table(fid, options, 'Phonon_Model');
  
  fprintf(fid, '<p>Here is some additional information about the atom/molecule physical properties\n');
  % display some information about the system (energy, structure, etc...)
  if isfield(Phonon_Model.UserData, 'properties')
    toadd = fieldnames(Phonon_Model.UserData.properties);
    fprintf(fid, '<table style="text-align: left; width: 50%%;" border="1" cellpadding="2" cellspacing="2">\n');
    for index=1:numel(toadd)
      this = Phonon_Model.UserData.properties.(toadd{index});
      if isnumeric(this) && ndims(this) <= 2
        if numel(this) > 50
          this1 = this(1:30);
          this2 = this((end-20):end);
          fprintf(fid, '  <tr><td><b>%s</b></td><td>%s ... %s</td></tr>\n', toadd{index}, mat2str(this1), mat2str(this2));
        else
          fprintf(fid, '  <tr><td><b>%s</b></td><td>%s</td></tr>\n', toadd{index}, mat2str(this));
        end
      elseif ischar(this)
        fprintf(fid, '  <tr><td><b>%s</b></td><td>%s</td></tr>\n', toadd{index}, this');
      end
    end
    fprintf(fid, '</table></p>\n');
  end
  fprintf(fid, '<p>End Date: %s.<br>\n', datestr(now));
  fprintf(fid, 'Time elapsed: %g [s]</p>\n', options.duration);
  if isfield(options, 'cite')
    fprintf(fid, '<p><b>Please cite:</b><br><pre>');
    fprintf(fid, '%s\n', options.cite{:});
    fprintf(fid, '</pre></p>\n');
  end
  
% ==============================================================================
function Phonon_DOS = sqw_phonons_htmlreport_dos(fid, options, object)
  % display the vDOS. Model must have been evaluated once to compute DOS
  Phonon_DOS = sqw_phonon_dos(object);
  [thermo, fig] = sqw_thermochemistry(object, 1:1000, 'plot');
  if isfield(object.UserData, 'DOS') && ~isempty(object.UserData.DOS)
    Phonon_DOS = object.UserData.DOS;
    if isfield(thermo, 'entropy'),          Thermo_S = thermo.entropy; end
    if isfield(thermo, 'internal_energy'),  Thermo_U = thermo.internal_energy; end
    if isfield(thermo, 'helmholtz_energy'), Thermo_F = thermo.helmholtz_energy; end
    if isfield(thermo, 'heat_capacity'),    Thermo_Cv= thermo.heat_capacity; end
    
    fprintf(fid, '<h3><a name="dos"></a>The vibrational density of states (vDOS)</h3>\n');
    fprintf(fid, '<p>The vibrational density of states (aka phonon spectrum) is defined as the velocity auto-correlation function (VACF) of the particles.\n<br>%s<br>\n', object.Name);
    
    if ~isempty(fig)
      plot2svg(fullfile(options.target,     'Phonon_ThermoChem.svg'), fig);
      print(fig,  fullfile(options.target,  'Phonon_ThermoChem.png'),'-dpng');
      print(fig,  fullfile(options.target,  'Phonon_ThermoChem.pdf'),'-dpdf');
      saveas(fig, fullfile(options.target,  'Phonon_ThermoChem.fig'),  'fig');
      close(fig);
      sqw_phonons_htmlreport_table(fid, options,'Phonon_ThermoChem');
    end
    
    properties = {'Phonon_DOS','Thermo_S','Thermo_U','Thermo_F','Thermo_Cv'};
    for index=1:numel(properties)
        name = properties{index};
        try
          val = eval(name);
          builtin('save', fullfile(options.target,  [ name '.mat']), name);
          save(val, fullfile(options.target, [ name '.dat']), 'dat data');
          save(val, fullfile(options.target, [ name '.h5']), 'mantid');
        catch
          disp([ mfilename ': ERROR when exporting ' name ' into HTML ' ])
          continue
        end
        switch name
        case 'Phonon_DOS'
          fprintf(fid, 'The phonon spectrum is available below:</p><br>\n');
        case 'Thermo_S'
          fprintf(fid, 'The entropy S=-dF/dT [J/K/mol] is available below:</p><br>\n');
        case 'Thermo_U'
          fprintf(fid, 'The internal energy U [J/mol] is available below:</p><br>\n');
        case 'Thermo_F'
          fprintf(fid, 'The Helmholtz free energy F=U-TS=-kT lnZ [J/mol] is available below:</p><br>\n');
        case 'Thermo_Cv'
          fprintf(fid, 'The molar specific heat at constant volume Cv=dU/dT [J/K/mol] is available below:</p><br>\n');
        end
        sqw_phonons_htmlreport_table(fid, options,[ name ]);
    
    end

  else Phonon_DOS=[];
  end

% ==============================================================================
function [maxFreq, object] = sqw_phonons_htmlreport_max_spectrum(fid, options, object)
  % compute the max energy of phonons, used for further evaluations
  qh=linspace(0.01,1.5,10);qk=qh; ql=qh; w=linspace(0.01,100,11);
  data=iData(object,object.p,qh,qk,ql',w);
  % search for a maxFreq item in the model.UserData
  if isfield(object.UserData, 'maxFreq')
    maxFreq = object.UserData.maxFreq*1.2;
  else
    maxFreq = max(w);
  end
  maxFreq = max(maxFreq(:));  % in case the max energies are given per mode
 
% ==============================================================================
function Phonon_kpath = sqw_phonons_htmlreport_kpath(fid, options, object, maxFreq)
  % generate dispersion along principal axes
  is_dos_there = isfield(object.UserData,'DOS') && ~isempty(object.UserData.DOS) ...
      && isa(object.UserData.DOS, 'iData') && prod(size(object.UserData.DOS)) > 10000;
  if ~is_dos_there, object.UserData.DOS = []; end
  [Phonon_kpath,kpath,fig] = sqw_kpath(object, 0, [0.01 maxFreq],'plot');
  if isfinite(max(Phonon_kpath)) && max(Phonon_kpath)
    Phonon_kpath = log10(Phonon_kpath/max(Phonon_kpath)); 
  else Phonon_kpath = log10(Phonon_kpath); end
  index=~isfinite(Phonon_kpath);
  Phonon_kpath(index) = nan;

  if ~isempty(Phonon_kpath)
    saveas(fig, fullfile(options.target, 'Phonon_kpath.fig'),'fig');
    plot2svg(fullfile(options.target, 'Phonon_kpath.svg'), fig);
    print(fig, '-dpdf',  fullfile(options.target, 'Phonon_kpath.pdf'));
    print(fig, '-dtiff', fullfile(options.target, 'Phonon_kpath.tiff'));
    close(fig);
    fprintf(fid, '<h3><a name="kpath"></a>The dispersion along principal directions</h3>\n');
    fprintf(fid, 'The dispersion curves along the principal axes is shown in log10 scale.<br>\n');
    builtin('save', fullfile(options.target, 'Phonon_kpath.mat'), 'Phonon_kpath');
    save(Phonon_kpath, fullfile(options.target, 'Phonon_kpath.png'), 'png data');
    save(Phonon_kpath, fullfile(options.target, 'Phonon_kpath.fits'), 'fits');
    save(Phonon_kpath, fullfile(options.target, 'Phonon_kpath.dat'), 'dat data');
    save(Phonon_kpath, fullfile(options.target, 'Phonon_kpath.h5'), 'mantid');
    save(Phonon_kpath, fullfile(options.target, 'Phonon_kpath.x3d'), 'x3d','tight auto');
    sqw_phonons_htmlreport_table(fid, options, 'Phonon_kpath');
  else close(fig); end
  
  if isfield(object.UserData,'FREQ') && isfield(object.UserData,'HKL')
    fprintf(fid, 'The HKL and Frequencies along the principal axes is available as:<br>\n');
    for f={'HKL','FREQ'}
      if isfield(object.UserData,f{1})
        s=object.UserData.(f{1}); 
        save(fullfile(options.target, [ 'Phonon_kpath_' f{1} '.dat' ]),  '-ascii', 's');
        fprintf(fid, '<a href="%s">%s</a><br>\n', ...
          fullfile(options.target, [ 'Phonon_kpath_' f{1} '.dat' ]), ...
          [ 'Phonon_kpath_' f{1} '.dat' ]);
      end
    end
    
  end

% ==============================================================================
function [Phonon_HKLE, object] = sqw_phonons_htmlreport_eval_4D(fid, options, object, maxFreq)
  % generate dispersion in 4D
  qh=linspace(0.01,1.5,30);qk=qh; ql=qh; w=linspace(0.01,maxFreq*1.1,101);
  Phonon_HKLE=iData(object,object.p,qh,qk,ql',w);

  % these are the biggest files
  fprintf(fid, '<h3><a name="grid4d"></a>The dispersion in 4D S(hkl,w)</h3>\n');
  fprintf(fid, '<p>This 4D Model has been evaluated on the grid (for rendering purposes):<br>\n');
  fprintf(fid, '<ul><li>qh=[%g:%g] with %i values (QH in rlu)</li>\n', ylim(Phonon_HKLE), size(Phonon_HKLE,1)); % axis1
  fprintf(fid, '<li>qk=[%g:%g] with %i values (QK in rlu)</li>\n', xlim(Phonon_HKLE), size(Phonon_HKLE,2));     % axis2
  fprintf(fid, '<li>ql=[%g:%g] with %i values (QL in rlu)</li>\n', zlim(Phonon_HKLE), size(Phonon_HKLE,3));
  fprintf(fid, '<li>w=[%g:%g] with %i values (energy in meV)</li></ul></p>\n', clim(Phonon_HKLE), size(Phonon_HKLE,4));
  fprintf(fid, '<p>The model parameters used for this evaluation are:<br>\n');
  fprintf(fid, '<pre>%s</pre><br></p>\n', class2str(' ', Phonon_HKLE.Parameters));
  
  fprintf(fid, 'Load the HKLE Data set under <a href="http://ifit.mccode.org">Matlab/iFit</a> with:<br>\n');
  fprintf(fid, '<ul><li>load(''<a href="%s">%s</a>'') <i>%% loads the 4D Phonon_HKLE iData object</i></li></ul><br>\n', 'Phonon_HKLE.mat', 'Phonon_HKLE.mat');
  
  builtin('save', fullfile(options.target, 'Phonon_HKLE.mat'), 'Phonon_HKLE');
  saveas(Phonon_HKLE, fullfile(options.target, 'Phonon_HKLE.dat'), 'dat data');
  saveas(Phonon_HKLE, fullfile(options.target, 'Phonon_HKLE.h5'), 'mantid');
  sqw_phonons_htmlreport_table(fid, options, 'Phonon_HKLE');

% ==============================================================================
function sqw_phonons_htmlreport_eval_3D(fid, options, object)
  % the S(0kl,w)
  x=object{1};
    
  Phonon_0KLE=object(1, :,:,:);
  clear object
  fprintf(fid, '<h3><a name="grid3d"></a>The dispersion in 3D S(h=0,k,l,w)</h3>\n');
  fprintf(fid, '<p>In order to view this 4D data set, we represent it on the QH~0 plane as a 3D volume data set. The intensity level is set as log10[S(QH=%g,QK,QL,w)].<br>\n', x(1));
    
  fprintf(fid, '<ul><li>qh=%g (QH in rlu)</li>\n', x(1)); % axis1
  fprintf(fid, '<li>qk=[%g:%g] with %i values (QK in rlu)</li>\n', xlim(Phonon_0KLE), size(Phonon_0KLE,2));     % axis2
  fprintf(fid, '<li>ql=[%g:%g] with %i values (QL in rlu)</li>\n', ylim(Phonon_0KLE), size(Phonon_0KLE,3));
  fprintf(fid, '<li>w=[%g:%g] with %i values (energy in meV, vertical)</li></ul></p>\n', zlim(Phonon_0KLE), size(Phonon_0KLE,4));
    
    fprintf(fid, '<p>The QH=%g data set is available in the folowing formats (log10 of the data set except for DAT and MAT files):<br>\n', x(1));
    
  builtin('save', fullfile(options.target, 'Phonon_0KLE.mat'), 'Phonon_0KLE');
  saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.dat'), 'dat data');
  Phonon_0KLE = log10(Phonon_0KLE);
  saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.png'), 'png', 'tight');
  saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.fig'), 'fig', 'plot3 tight');
  saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.vtk'));
  saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.mrc'));
  
  % the PDF export may crash in deployed version
  if ~isdeployed
    saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.pdf'), 'pdf', 'tight');
  end
  saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.h5'), 'mantid');
  
  % modify aspect ratio to fit in a cube for X3D/XHTML
  saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.xhtml'), 'xhtml','axes auto');
  saveas(Phonon_0KLE, fullfile(options.target, 'Phonon_0KLE.x3d'), 'x3d','axes auto');
  
  sqw_phonons_htmlreport_table(fid, options, 'Phonon_0KLE');

% ==============================================================================
function Phonon_powder = sqw_phonons_htmlreport_eval_powder(fid, options, object, maxFreq)
  % the S(hklw) radial average
  Phonon_powder = sqw_powder(object, [], [], [0 maxFreq]); 
  log_Phonon_powder=log(Phonon_powder);
  
  fprintf(fid, '<h3><a name="powder"></a>The powder average S(q,w)</h3>\n');
  fprintf(fid, 'The powder average S(q,w) is shown below:<br>\n');
    
  builtin('save', fullfile(options.target, 'Phonon_powder.mat'), 'Phonon_powder');
  saveas(log_Phonon_powder, fullfile(options.target, 'Phonon_powder.png'),'png data');
  saveas(log_Phonon_powder, fullfile(options.target, 'Phonon_powder.fits'),'fits');
  saveas(log_Phonon_powder, fullfile(options.target, 'Phonon_powder.tiff'),'tiff');
  saveas(log_Phonon_powder, fullfile(options.target, 'Phonon_powder.fig'), 'fig', 'tight');
  saveas(Phonon_powder, fullfile(options.target, 'Phonon_powder.dat'), 'dat data');
  saveas(log_Phonon_powder, fullfile(options.target, 'Phonon_powder.pdf'), 'pdf', 'tight view2');
  saveas(log_Phonon_powder, fullfile(options.target, 'Phonon_powder.svg'), 'svg', 'view2 tight');
  saveas(log_Phonon_powder, fullfile(options.target, 'Phonon_powder.x3d'), 'x3d', 'tight auto');
  saveas(Phonon_powder, fullfile(options.target, 'Phonon_powder.h5'), 'mantid');
  
  sqw_phonons_htmlreport_table(fid, options, 'Phonon_powder');

function sqw_phonons_htmlreport_download(fid, options)
  fprintf(fid, '<h2><a name="zip">Download</h2>\n');
  fprintf(fid, 'You can download the whole content of this report (data, plots, logs, ...) from<br>\n');
  [~,short_path] = fileparts(options.target);
  fprintf(fid, '<ul><li><a href="%s">%s</a></li></ul>', fullfile('..',[ short_path '.zip' ]), [ short_path '.zip' ]);
  fprintf(fid, 'You should then open the "sqw_phonons.html" file therein with a browser (Firefox, Chrome...)<br>\n');
  
  % close HTML document
  fprintf(fid, '<hr>Date: %s.<br>\n', datestr(now));
  fprintf(fid, 'Powered by <a href="http://ifit.mccode.org">iFit</a> E. Farhi (c) 2016.\n</body></html>');
  
  % create a simple README file
  freadme = fopen(fullfile(options.target,'README.txt'),'w');
  fprintf(freadme, 'Open the sqw_phonons.html file in this directory. It contains all you need. The Phonon Model (iFunc) is stored in the Phonon_Model.mat Matlab binary file.\n');
  fclose(freadme);
  
  % create ZIP of document
  zip(options.target, options.target);
  
  if isdeployed || ~usejava('jvm') || ~usejava('desktop')
    disp([ 'sqw_phonons: HTML report created as ' options.report ]);
  else
    disp([ 'sqw_phonons: HTML report created as <a href="' options.report '">' options.report '</a>' ]);
  end



function sqw_phonons_htmlreport(filename, step, options, data, object)

if isempty(filename), return; end
if ~options.htmlreport, return; end

if strcmp(step, 'init')
  fid = fopen(filename, 'w');
else
  fid = fopen(filename, 'a');
end

if fid < 0, return; end

switch step
case 'init'
  % retrieve a few external resources
  copyfile(fullfile(ifitpath,'Docs','images','iFit-logo.png'), options.target);
  % http://ifit.mccode.org
  copyfile(fullfile(ifitpath,'Docs','images','ase256.png'), options.target);
  % https://wiki.fysik.dtu.dk/ase
  switch options.calculator
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
  case 'QUANTUMESPRESSO'
  logo='logo_qe.jpg';
  link='http://www.quantum-espresso.org/';
  copyfile(fullfile(ifitpath,'Docs','images',logo), options.target);
  otherwise
  logo=''; link='';
  end
  % open the report HTML page and write header, title, date, ...
  fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n');
  fprintf(fid, '<html>\n<head>\n<title>%s: %s with %s</title>\n</head>\n', ...
    'sqw_phonons', options.configuration, options.calculator);
  fprintf(fid, '<body><div style="text-align: center;">\n');
  fprintf(fid, '<a href="http://ifit.mccode.org"><img title="ifit.mccode.org" src="iFit-logo.png" align="middle" height=100></a>\n');
  fprintf(fid, '<a href="https://wiki.fysik.dtu.dk/ase" title="wiki.fysik.dtu.dk/ase"><img src="ase256.png" align="middle" height=100></a></br>\n');
  fprintf(fid, '<h1>%s: %s with %s</h1></div>\n', 'sqw_phonons', options.configuration, options.calculator);
  % general introduction
  fprintf(fid, 'Start Date: %s.<br>\nComputer: %s.<br>\n', datestr(now), computer);
  fprintf(fid, 'Stored in: <a href="%s">%s</a><br>\n', options.target, options.target);
  fprintf(fid, '<p>This page presents the results of the estimate of phonon dispersions in a single crystal, using a DFT code (the "calculator") in the background. From the initial atomic configuration (geometry), each atom in the lattice cell is displaced by a small quantity. The displaced atom then sustains a, so called Hellmann-Feynman, restoring force to come back to the stable structure. The dynamical matrix is obtained from these forces, and its eigen-values are the energies of the vibrational modes in the crystal.</p>\n');
  fprintf(fid, '<p>This computational resource is provided by <a href="http://ifit.mccode.org">iFit</a>, with the <b>sqw_phonon</b> Model, which itself makes use of the <a href="https://wiki.fysik.dtu.dk/ase">Atomistic Simulation Environment (ASE)</a>.\n');
  fprintf(fid, '<p>This report summarizes the initial crystal geometry and the calculator configuration. When the simulation ends succesfully, the lower part presents the dispersion curves as plots and data sets, the model used (as a Matlab object), and the density of states. These results correspond to the coherent inelastic vibrational modes.</p>\n');

  % append system configuration
  fprintf(fid, '<h2>Atom/molecule configuration</h2>\n');
  fprintf(fid, '<p>In this section, we present the crystallographic lattice cell used for the computation, both as plots and structure files for use with other software.</p>\n');
  fprintf(fid, 'Initial description: %s<br>\n', options.configuration);
  fprintf(fid, '<p>\n');
  % append links to configuration files and images
  for index={ 'configuration.html You can rotate the model (left mouse button), zoom (right mouse button), and translate (middle mouse button).', ...
              'configuration.png', ...
              'configuration.cif Crystallographic Information File (<a href="https://en.wikipedia.org/wiki/Crystallographic_Information_File">CIF</a>) which you can view with <a href="http://jmol.sourceforge.net/">JMol</a>, <a href="http://www.ks.uiuc.edu/Research/vmd/"VMD</a>, <a href="http://www.cgl.ucsf.edu/chimera/">Chimera</a>, <a href="http://rasmol.org/">RasMol</a>, ...', ...
              'configuration.pdb Protein Data Bank (<a href="http://www.rcsb.org/pdb/">PDB</a>) which you can view with <a href="http://jmol.sourceforge.net/">JMol</a>, <a href="http://www.ks.uiuc.edu/Research/vmd/"VMD</a>, <a href="http://www.cgl.ucsf.edu/chimera/">Chimera</a>, <a href="http://rasmol.org/">RasMol</a>, ...', ...
              'configuration.etsf <a href="http://www.etsf.eu/fileformats">European Theoretical Spectroscopy Facility</a> file format for DFT codes', ...
              'configuration_SHELX.res <a href="http://shelx.uni-ac.gwdg.de/SHELX/">ShelX</a> file format', ...
              'configuration_VASP POSCAR geometry file for <a href="https://www.vasp.at/">VASP</a>', ...
              'configuration.eps Encapsulated postscript', ...
              'configuration.pov Scene for <a href="http://www.povray.org/">Pov-Ray</a>', ...
              'configuration.x3d Scene for <a href="http://castle-engine.sourceforge.net/view3dscene.php">view3dscene</a>, <a href="http://www.instantreality.org/">InstantPlayer</a>, <a href="http://freewrl.sourceforge.net/">FreeWRL</a>'  }
    [index1, index2] = strtok(index{1});
    if strcmp(index1, 'configuration.png')
      fprintf(fid, '<div style="text-align:center"><img src="%s" align="middle" title="%s"></div><br>\n', index1, index1);
    elseif strcmp(index1, 'configuration.html')
      fprintf(fid, '<div style="text-align:center"><iframe src="%s" onload="this.style.height=this.contentDocument.body.scrollHeight +''px'';" align="middle"></iframe><br>%s (<a href="%s" target=_blank>open in external window</a>)<br></div><br>\n', index1, index2, index1);
    else
      fprintf(fid, '[ <a href="%s">%s</a> ] %s<br>\n', index1, index1, index2);
    end

  end
  fprintf(fid, '</p>\n');
  % append calculator configuration
  fprintf(fid, '<h2>Calculator configuration</h2>\n');
  % general introduction
  fprintf(fid, 'This is a quick desription of the calculator used to compute the lattice energy and forces.<br>\n');
  if ~isempty(logo)
    fprintf(fid, 'We are using the calculator: <a href="%s"><img src="%s" height="80" align="middle"><b>%s</b></a> (<a href="%s">%s</a>).<br>\n', ...
      link, logo, upper(options.calculator), link, link);
  else
    fprintf(fid, 'The <b>%s</b> calculator call is:<br>\n', upper(options.calculator));
  end
  fprintf(fid, '  <ul><li><code>%s</code></li></ul>\n', data); % 'calc'
  % clean 'options' from empty and 0 members
  op           = options;
  op.htmlreport= 0;
  for f=fieldnames(op)'
    if isempty(op.(f{1})) || (isscalar(op.(f{1})) && op.(f{1}) == 0)
      op = rmfield(op, f{1});
    end
  end
  fprintf(fid, 'The calculator configuration:<br><div style="margin-left: 40px;"><pre>%s</pre></div>\n', class2str(' ',op));
  % indicate that we are computing
  fprintf(fid, '<b>Computing... (be patient)</b>\n');
  
case 'done'
  % indicate evaluated model, and print grid used
  fprintf(fid, '<h2>Computation completed</h2>\n');
  fprintf(fid, '<b>WELL DONE</b><br><p>The computation has performed correctly. The forces and dynamical matrix of the lattice vibrations has been determined.</p><br>\n');
  fprintf(fid, 'End Date: %s.<br>\n', datestr(now));
  fprintf(fid, 'Time elapsed: %g [s]<br>\n', options.duration);
  fprintf(fid, 'Please cite:<br><pre>');
  fprintf(fid, '%s\n', data{:});
  fprintf(fid, '</pre>\n');
  
case 'plot'
  % present the final object
  fprintf(fid, '<h2>Results</h2>\n');
  fprintf(fid, '<h3>The Phonon dispersion Model</h3>\n');
  Phonons = object;
  builtin('save', fullfile(options.target, 'iFunc_Phonons.mat'), 'Phonons');
  fprintf(fid, '<p>The results are stored into a 4D <a href="http://ifit.mccode.org?iFunc.html">iFunc</a> object containing the dynamical matrix. This is a Matlab workspace (MAT-file).');
  fprintf(fid, '<ul><li><a href="%s">%s</a></li></ul>\n', ...
    'iFunc_Phonons.mat', 'iFunc_Phonons.mat');
  fprintf(fid, 'Load the Model under <a href="http://ifit.mccode.org">Matlab/iFit</a> with (this also works with the <a href="http://ifit.mccode.org/Install.html">standalone version of iFit</a> which does <b>not</b> require any Matlab license and installs on most systems): <ul><li>load(''<a href="%s">%s</a>'') <i>%% creates a "Phonons" iFunc Model</i></li></ul>\n', 'iFunc_Phonons.mat', 'iFunc_Phonons.mat');
  fprintf(fid, 'Define axes for the evaluation grid in 4D, for instance:<ul><li>qh=linspace(0.01,.5,50); qk=qh; ql=qh; w=linspace(0.01,50,51);</li></ul>\n');
  fprintf(fid, 'Evaluate the Phonons under Matlab/iFit with: <ul><li>iData(Phonons, [], qh, qk, ql, w) <i>%% evaluates the "Phonons" onto the grid, with default parameters, and return an iData object</i></li></ul></p>\n');
  
  % append evaluated plots and link to data sets (VTK, MCR, images, ...)
  fprintf(fid, '<h3>Evaluating the Phonon dispersion Model onto a grid</h3>\n');
  if isempty(data)
    qh=linspace(0.01,.5,20);qk=qh; ql=qh; w=linspace(0.01,100,51);
    data=iData(object,[],qh,qk,ql,w);
  end
  Phonons_HKLE = data;
  builtin('save', fullfile(options.target, 'iData_Phonons.mat'), 'Phonons_HKLE');
  saveas(Phonons_HKLE, fullfile(options.target, 'Phonons_HKLE.dat'));
  saveas(Phonons_HKLE, fullfile(options.target, 'Phonons_HKLE.h5'), 'mantid');
  clear Phonons_HKLE
  fprintf(fid, '<p>This 4D Model has been evaluated on the grid (for rendering purposes):<br>\n');
  fprintf(fid, '<ul><li>qh=[%g:%g] with %i values (QH in rlu)</li>\n', ylim(data), size(data,1)); % axis1
  fprintf(fid, '<li>qk=[%g:%g] with %i values (QK in rlu)</li>\n', xlim(data), size(data,2));     % axis2
  fprintf(fid, '<li>ql=[%g:%g] with %i values (QL in rlu)</li>\n', zlim(data), size(data,3));
  fprintf(fid, '<li>w=[%g:%g] with %i values (energy in meV)</li></ul>\n', clim(data), size(data,4));
  
  fprintf(fid, 'Load the HKLE Data set under <a href="http://ifit.mccode.org">Matlab/iFit</a> with:<br>\n');
  fprintf(fid, '<ul><li>load(''<a href="%s">%s</a>'') <i>%% loads the 4D Phonons_HKLE iData object</i></li>\n', 'iData_Phonons.mat', 'iData_Phonons.mat');
  fprintf(fid, '<li>[ <a href="%s">%s</a> ] is a flat text file which contains axes and the 4D data set. You will have to reshape the matrix after reading the contents.</li>\n', 'Phonons_HKLE.dat', 'Phonons_HKLE.dat');
  fprintf(fid, '<li>[ <a href="%s">%s</a> ] a NeXus/HDF5 data file to be opened with e.g. <a href="http://www.mantidproject.org/Main_Page">Mantid</a> or <a href="http://www.hdfgroup.org/hdf-java-html/hdfview">hdfview</a>.</li></ul></p>\n', 'Phonons_HKLE.h5', 'Phonons_HKLE.h5');
  fprintf(fid, '<p>In order to view this 4D data set, we represent it on the QH~0 plane as a 3D volume data set.<br>\n');
  x=data{1};
  fprintf(fid, '<ul><li>qh=%g (QH in rlu)</li>\n', x(1)); % axis1
  fprintf(fid, '<li>qk=[%g:%g] with %i values (QK in rlu)</li>\n', xlim(data), size(data,2));     % axis2
  fprintf(fid, '<li>ql=[%g:%g] with %i values (QL in rlu)</li>\n', zlim(data), size(data,3));
  fprintf(fid, '<li>w=[%g:%g] with %i values (energy in meV, vertical)</li></ul></p>\n', clim(data), size(data,4));
  % TODO: add details on axes and isosurface
  % TODO: add DOS
  data=-log(data(1, :,:,:));
  saveas(data, fullfile(options.target, 'Phonons.png'));
  saveas(data, fullfile(options.target, 'Phonons.vtk'));
  saveas(data, fullfile(options.target, 'Phonons.mrc'));
  saveas(data, fullfile(options.target, 'Phonons3D.dat'));
  saveas(data, fullfile(options.target, 'Phonons3D.h5'), 'mantid');
  % modify aspect ratio to fit in a cube
  data{1}=linspace(0,1,size(data,1));
  data{2}=linspace(0,1,size(data,2));
  data{3}=linspace(0,1,size(data,3));
  
  fprintf(fid, '<p>The QH=0 data set is available in the folowing formats:<br>\n');
  fprintf(fid, '<ul><li>[ <a href="%s">%s</a> ] Visualization Toolkit (VTK) file which can be viewed with <a href="http://www.paraview.org/">ParaView</a>, <a href="http://code.enthought.com/projects/mayavi/">Mayavi2</a>, <a href="https://wci.llnl.gov/simulation/computer-codes/visit/executables">VisIt</a> and <a href="https://www.slicer.org/">Slicer4</a>.</li>\n', 'Phonons.vtk', 'Phonons.vtk');
  
  fprintf(fid, '<li>[ <a href="%s">%s</a> ] electron density map format (MRC) which can be viewed with PyMol, <a href="http://www.ks.uiuc.edu/Research/vmd/">VMD</a>, <a href="http://www.cgl.ucsf.edu/chimera/">Chimera</a>, <a href="http://www.yasara.org/">Yasara</a>.</li>\n', 'Phonons.mrc', 'Phonons.mrc');
  fprintf(fid, '<li>[ <a href="%s">%s</a> ] a flat text file which contains axes and the 3D data set. You will have to reshape the matrix after reading the contents.</li>\n', 'Phonons3D.dat', 'Phonons3D.dat');
  fprintf(fid, '<li>[ <a href="%s">%s</a> ] a NeXus/HDF5 data file to be opened with e.g. <a href="http://www.mantidproject.org/Main_Page">Mantid</a> or <a href="http://www.hdfgroup.org/hdf-java-html/hdfview">hdfview</a>.</li></ul>\n', 'Phonons3D.h5', 'Phonons3D.h5');
  fprintf(fid, '</p>\n');
  
  fprintf(fid, 'Now, here are a few representations of the 3D data set<br>\n');
  saveas(data, fullfile(options.target, 'Phonons.xhtml'));
  
  fprintf(fid, '<img src="%s" title="%s" align="middle"><br>\n', 'Phonons.png', 'Phonons.png');
  fprintf(fid, '<iframe src="%s" align="middle" width="700" height="1000"></iframe><br>The blue axis is the Energy (meV), the red axis is the QK (rlu), the green axis is QL (rlu).<br>\nYou can rotate the model (left mouse button), zoom (right mouse button), and translate (middle mouse button). <br>(<a href="%s" target=_blank>open in external window</a>)<br>\n', 'Phonons.xhtml', 'Phonons.xhtml');
  
case 'error'
  fprintf(fid, '<h2>ERROR: %s</h2>\n', data);
  fprintf(fid, [ mfilename ' ' options.configuration ' ' options.calculator ' <b>FAILED</b><br>' ]);
  sqw_phonons_htmlreport(filename, 'end');
case 'end'
  % close HTML document
  fprintf(fid, '<hr>Date: %s.<br>\n', datestr(now));
  fprintf(fid, 'Powdered by <a href="http://ifit.mccode.org>iFit</a> E. Farhi (c) 2016.\n');
  if isdeployed || ~usejava('jvm') || ~usejava('desktop')
    disp([ 'HTML report created in ' options.target ]);
  else
    disp([ 'HTML report created in <a href="' fullfile(options.target,'index.html') '">' options.target '</a>' ]);
  end
end

fclose(fid);


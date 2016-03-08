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
  copyfile(fullfile(ifitpath,'Models','images','iFit-logo.png'), options.target);
  % http://ifit.mccode.org
  copyfile(fullfile(ifitpath,'Models','images','ase256.png'), options.target);
  % https://wiki.fysik.dtu.dk/ase
  switch options.calculator
  case 'GPAW'
  logo='logo-gpaw.png';
  link='http://wiki.fysik.dtu.dk/gpaw';
  copyfile(fullfile(ifitpath,'Models','images',logo), options.target);
  case 'NWCHEM'
  logo='nwchem.png';
  link='http://www.nwchem-sw.org/';
  copyfile(fullfile(ifitpath,'Models','images',logo), options.target);
  case 'ELK'
  logo='elk.png';
  link='http://elk.sourceforge.net';
  copyfile(fullfile(ifitpath,'Models','images',logo), options.target);
  case {'DACAPO','JACAPO'}
  logo='jacapo.png';
  link='http://wiki.fysik.dtu.dk/dacapo';
  copyfile(fullfile(ifitpath,'Models','images',logo), options.target);
  case 'ABINIT'
  logo='abinit.png';
  link='http://www.abinit.org/';
  copyfile(fullfile(ifitpath,'Models','images',logo), options.target);
  case 'EMT'
  logo='emt.png';
  link='https://wiki.fysik.dtu.dk/ase/ase/calculators/emt.html';
  copyfile(fullfile(ifitpath,'Models','images',logo), options.target);
  case 'QUANTUMESPRESSO'
  logo='logo_qe.jpg';
  link='http://www.quantum-espresso.org/';
  copyfile(fullfile(ifitpath,'Models','images',logo), options.target);
  otherwise
  logo=''; link='';
  end
  % open the report HTML page and write header, title, date, ...
  fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n');
  fprintf(fid, '<html>\n<head>\n<title>%s: %s with %s</title>\n</head>\n', ...
    'sqw_phonons', options.configuration, options.calculator);
  fprintf(fid, '<body>\n');
  fprintf(fid, '<a href="http://ifit.mccode.org"><img src="iFit-logo.png" align="middle" height=100></a>\n');
  fprintf(fid, '<a href="https://wiki.fysik.dtu.dk/ase"><img src="ase256.png" align="middle" height=100></a></br>\n');
  fprintf(fid, '<h1>%s: %s with %s</h1>\n', 'sqw_phonons', options.configuration, options.calculator);
  fprintf(fid, 'Start Date: %s.<br>\nComputer: %s.<br>\n', datestr(now), computer);
  fprintf(fid, 'Stored in: <a href="%s">%s</a><br>\n', options.target, options.target);

  % append system configuration
  fprintf(fid, '<h2>Atom/molecule configuration</h2>\n');
  fprintf(fid, 'Initial description: %s<br>\n', options.configuration);
  fprintf(fid, '<p>\n');
  % append links to configuration files and images
  for index={ 'configuration.html','configuration.png', ...
              'configuration.cif', 'configuration.pdb', ...
              'configuration.etsf','configuration_SHELX.res', 'configuration_VASP', ...
              'configuration.eps', 'configuration.pov', 'configuration.x3d'  }
    this = fullfile(options.target, index{1});
    if ~isempty(dir(this))
      if strcmp(index{1}, 'configuration.png')
        fprintf(fid, '<img src="%s" align="middle"><br>\n', index{1}, index{1}, index{1});
      elseif strcmp(index{1}, 'configuration.html')
        fprintf(fid, '<iframe src="%s" onload="this.style.height=this.contentDocument.body.scrollHeight +''px'';" align="middle"></iframe><br>\n', index{1});
      else
        fprintf(fid, '[ <a href="%s">%s</a> ]<br>\n', index{1}, index{1});
      end
    end
  end
  fprintf(fid, '</p>\n');
  % append calculator configuration
  fprintf(fid, '<h2>Calculator configuration</h2>\n');
  if ~isempty(logo)
    fprintf(fid, 'We are using <a href="%s"><img src="%s" height="80"><b>%s</b></a> (<a href="%s">%s</a>).<br>\n', ...
      link, logo, upper(options.calculator), link, link);
  else
    fprintf(fid, 'We are using <b>%s</b>.<br>\n', upper(options.calculator));
  end
  fprintf(fid, '  %s<br>\n', data); % 'calc'
  % clean 'options' from empty and 0 members
  op = options;
  for f=fieldnames(op)
    if isempty(op.(f{1})) || op.(f{1}) == 0
      op = rmfield(op, f{1});
    end
  end
  fprintf(fid, 'Calculator configuration:<br><p><pre>%s</pre></p>\n', class2str(op));
  % indicate that we are computing
  fprintf(fid, '<b>Computing... (be patient)</b>\n');
  
case 'done'
  % indicate evaluated model, and print grid used
  fprintf(fid, '<h2>Computation completed</h2>\n');
  fprintf(fid, 'End Date: %s.<br>\n', datestr(now));
  fprintf(fid, 'Time elapsed: %g [s]<br>\n', options.duration);
  fprintf(fid, 'Please cite:<br><pre>');
  fprintf(fid, '%s\n', data{:});
  fprintf(fid, '</pre>\n');
  
case 'plot'
  % present the final object
  fprintf(fid, '<h2>Results</h2>\n');
  Phonons = object;
  builtin('save', fullfile(options.target, 'iFunc_Phonons.mat'), Phonons);
  fprintf(fid, '<p>The results are stored into an <a href="http://ifit.mccode.org?iFunc.html">iFunc</a> object containing the dynamical matrix.');
  fprintf(fid, '<ul><li><a href="%s">%s</a></li></ul>\n', ...
    'iFunc_Phonons.mat', 'iFunc_Phonons.mat');
  fprintf(fid, 'Load the Model under <a href="http://ifit.mccode.org">Matlab/iFit</a> with: <ul><li>Phonons=load(''<a href="%s">%s</a>'')</li></ul>\n', 'iFunc_Phonons.mat', 'iFunc_Phonons.mat');
  fprintf(fid, 'Define axes for the evaluation grid in 4D, for instance:<ul><li>qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);</li></ul>\n');
  fprintf(fid, 'Evaluate the Phonons under Matlab/iFit with: <ul><li>iData(Phonons, [], qh, qk, ql, w)</li></ul></p>\n');
  
  % append evaluated plots and link to data sets (VTK, MCR, images, ...)
  if isempty(data)
    qh=linspace(0.01,.5,10);qk=qh; ql=qh; w=linspace(0.01,150,11);
    data=iData(object,[],qh,qk,ql,w);
  end
  Phonons_HKLE = data; clear data;
  fprintf(fid, '<p>This Model has been evaluated on grid:<br>\n');
  fprintf(fid, '<ul><li>QH=[%g:%g]</li>\n', ylim(data)); % axis1
  fprintf(fid, '<li>QK=[%g:%g]</li>\n', xlim(data));     % axis2
  fprintf(fid, '<li>QL=[%g:%g]</li></ul>\n', zlim(data));
  builtin('save', fullfile(options.target, 'iData_Phonons.mat'), Phonons_HKLE);
  clear Phonons_HKLE
  fprintf(fid, 'Load the HKLE Data set under <a href="http://ifit.mccode.org">Matlab/iFit</a> with:<br>\n');
  fprintf(fid, '<ul><li>Phonons_HKLE=load(''<a href="%s">%s</a>'')</li></ul></p>\n', 'iData_Phonons.mat', 'iData_Phonons.mat');
  fprintf(fid, '<p>In order to view this 4D data set, we represent it on the QH=0 plane.<br>\n');
  data=-log(data(1, :,:,:));
  saveas(data, fullfile(options.target, 'Phonons.png'));
  saveas(data, fullfile(options.target, 'Phonons.pdf'), '', 'plot3 view3');
  saveas(data, fullfile(options.target, 'Phonons.xhtml'));
  saveas(data, fullfile(options.target, 'Phonons.vtk'));
  fprintf(fid, '<img src="%s" align="middle"><br>\n', 'Phonons.png');
  fprintf(fid, '<iframe src="%s" onload="this.style.height=this.contentDocument.body.scrollHeight +''px'';" align="middle"></iframe><br>\n', 'Phonons.xhtml');
  fprintf(fid, '[ <a href="%s">%s</a> ]<br>\n', 'Phonons.pdf', 'Phonons.pdf');
  fprintf(fid, '[ <a href="%s">%s</a> ]<br>\n', 'Phonons.vtk', 'Phonons.vtk');
  fprintf(fid, '</p>\n');
case 'error'
  fprintf(fid, '<h2>ERROR: %s</h2>\n', data);
  fprintf(fid, [ mfilename ' ' options.configuration ' ' options.calculator ' FAILED<br>' ]);
  sqw_phonons_htmlreport(filename, 'end');
case 'end'
  % close HTML document
  fprintf(fid, 'End Date: %s.<br>\n', datestr(now));
  fprintf(fid, '<hr>\nPowdered by <a href="http://ifit.mccode.org>iFit</a> E. farhi (c) 2016.\n');
  fclose(fid)
end

fclose(fid);


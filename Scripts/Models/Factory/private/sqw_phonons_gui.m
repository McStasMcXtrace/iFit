function options = sqw_phonons_gui(configuration, options, status)

  % pop-up a simple dialog requesting for:
  %  * configuration
  %  * calculator
  %  * metal, insulator, semiconductor
  %  * supercell
  %  * kpoints
  % and will select autoplot, 'dos'
  doc(iData,'Models.html#mozTocId990577');  % doc for Phonons
  calcs = 'EMT';
  for index={'gpaw','elk','jacapo','nwchem','abinit','quantumespresso','quantumespresso_ase','vasp'};
    if ~isempty(status.(index{1})), calcs = [ calcs ', ' upper(index{1}) ]; end
  end
  calcs = strrep(calcs, '_','\_');
  NL = sprintf('\n');
  prompt = { [ '{\bf Atom/molecule/system configuration}' NL 'a CIF/PDB/POSCAR/... name or e.g. bulk("Cu", "fcc", a=3.6, cubic=True),' NL  'molecule("H2O"), or nanotube(6, 0, length=4). ' NL 'Documentation at {\color{blue}http://ifit.mccode.org/Models.html}' ], ...
  [ '{\bf Calculator}' NL 'One of ' calcs NL '{\color{red}BEWARE the computation may be LONG (days)}. We recommend QuantumEspresso and ABINIT.' ], ...
  [ '{\bf Smearing}' NL 'metal, semiconductor, insulator or left empty. Use e.g. 0.3 eV for conductors, or a small value such as 0.01 to help SCF convergence. You may use "auto" with Elk. ' ], ...
  [ '{\bf Cut-off energy for wave-functions}' NL 'Leave as 0 for the default, or specify a cut-off in eV, e.g. 500 eV for fast estimate, 1000 for ABINIT, 1500 or 2000 eV for accurate results.' ], ...
  [ '{\bf K-Points}' NL 'Monkhorst-Pack grid which determines the K sampling (3-vector). 4 is the minimum for accurate computations, 6 is best. Use 1 or 2 for testing only (faster).' ], ...
   [ '{\bf Supercell}' NL 'The size of the repeated model = system*supercell (3-vector). Should be larger than k-points.' ], ...
   [ '{\bf Other options}' NL 'Such as mpi, nbands, nsteps, xc (default PBE), toldfe, raw' NL 'example: "mpi=4; nsteps=100"' NL 'Documentation at {\color{blue}http://ifit.mccode.org/Models.html}' ] };
  dlg_title = 'iFit: Model: Sqw phonons';
  defAns    = {configuration, options.calculator, options.occupations, num2str(options.ecut), ...
    num2str(options.kpoints), num2str(options.supercell), ''};
  num_lines = [ 1 1 1 1 1 1 1 ]';
  op.Resize      = 'on';
  op.WindowStyle = 'normal';   
  op.Interpreter = 'tex';
  answer = inputdlg(prompt, dlg_title, num_lines, defAns, op);
  if isempty(answer), 
    options = [];
    return; 
  end
  % extract results
  configuration      = answer{1};
  options.calculator = answer{2};
  options.occupations= answer{3};
  if ~isnan(str2double(options.occupations))
    options.occupations = str2double(options.occupations);
  end
  options.ecut       = str2num(answer{4});
  
  options.kpoints    = str2num(answer{5});
  options.supercell  = str2num(answer{6});
  % other options transfered to 'options'
  others             = str2struct(answer{7});
  if ~isempty(others)
    for f=fieldnames(others)
      if ~isfield(options, f{1}) || isempty(options.(f{1}))
        options.(f{1}) = others.(f{1});
      end
    end
  end
  options.autoplot   = 1;
  options.gui        = true;
  if strcmpi(options.calculator,'qe_ase') options.calculator='quantumespresso_ase'; 
  elseif strcmpi(options.calculator,'qe')         options.calculator='quantumespresso'; 
  end

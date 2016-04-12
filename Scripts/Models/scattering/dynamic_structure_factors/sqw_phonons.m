function signal=sqw_phonons(configuration, varargin)
% model=sqw_phonons(configuration, ..., options)
%
%   iFunc/sqw_phonons: computes phonon dispersions using the ASE.
%   A model which computes phonon dispersions from the forces acting between
%     atoms. The input argument is any configuration file describing the
%     material, e.g. CIF, PDB, POSCAR, ... supported by ASE.
%   This models can compute the coherent inelastic phonons dispersions for
%     any crystalline (powder or single crystal) material.
%   The phonon spectra is computed using one of the calculator supported by the
%   Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase>.
%   Supported calculators are:
%     EMT       Effective Medium Theory calculator (Al,Cu,Ag,Au,Ni,Pd,Pt,H,C,N,O)
%     GPAW      Real-space/plane-wave/LCAO PAW code
%     NWChem    Gaussian based electronic structure code
%     Dacapo    Plane-wave ultra-soft pseudopotential code
%     ELK       Full Potential LAPW code
%     ABINIT    Plane-wave pseudopotential code
%     QuantumEspresso Plane-wave pseudopotential code
%     VASP      Plane-wave PAW code (when installed and licensed)
%
%   We recommend QuantumEspresso, GPAW, ABINIT and Elk.
%   Refer to https://molmod.ugent.be/deltacodesdft for codes comparison.
%
%   The simplest usage is to call: sqw_phonons('gui') which displays an entry dialog box
%   and proceeds with a fully automatic computation, and plots final results.
%     
%   The calculators can be specified by just giving their name as a parameter, 
%   or using e.g. options.calculator='GPAW'. Except for EMT, other calculators must
%   be installed separately, and optionally specific pseudo-potentials. 
%     See https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html
%
%   When performing a model evaluation, the DOS is also computed and stored
%     when the options 'dos' is specified during the creation. The DOS is only 
%     computed during the first evaluation, and is stored in model.UserData.DOS 
%     as an iData object. Subsequent evaluations are faster.
%
%   Benchmarks indicate that, for phonon dispersions:
%   * QuantumEspresso is the fastest, with excellent parallelization and accuracy.
%   * ABINIT, VASP are also fast codes.
%   * The all-electrons Elk code is about 10 times slower than QuantumEspresso.
%   * Using k-points=3 is both fast and ensures a reasonable accuracy in phonons
%   * Phonon dispersions are not too sensitive on the energy cut-off. 340 eV is good.
%   * Elk and ABINIT can not handle large supercells without recompiling.
%   * NWChem and Elk are not sensitive to the energy cut-off.
%
% The arguments for the model creation should be:
%
% input (model creation):
% configuration: file name to an existing material configuration
%   Any A.S.E supported format can be used (POSCAR, CIF, SHELX, PDB, ...). 
%     See <https://wiki.fysik.dtu.dk/ase/ase/io.html#module-ase.io>
%   Alternatively, the 'bulk','molecule', 'crystal', and 'nanotube' ASE constructors can be
%   used, using the Python syntax, e.g. 
%       'bulk("Si", "diamond", a=5.4)'
%       'bulk("Cu", "fcc", a=3.6, cubic=True)'
%       'molecule("H2O")'
%       'nanotube(6, 0, length=4)'
%       'crystal(["Na", "Cl"], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
%          cellpar=[5.64, 5.64, 5.64, 90, 90, 90])'
%     See <https://wiki.fysik.dtu.dk/ase/ase/structure.html>
%         <https://wiki.fysik.dtu.dk/ase/tutorials/spacegroup/spacegroup.html>
%
% 'metal','insulator','semiconductor': indicates the type of occupation for 
%    electronic states, which sets smearing.
%
% 'dos': will compute the vibrational density of states.
% 'html' or 'report':   generate a full report of the simulation and export data.
% 'plot' or 'autoplot': plot the results after building the model.
%
% options: an optional structure with optional settings
%
% General options
%   options.target =path                   Where to store all files and FORCES
%     a temporary directory is created when not specified.
%   options.supercell=scalar or [nx ny nz] supercell size. Default is 3.
%   options.calculator=string              EMT,GPAW,Elk,NWChem,Dacapo,ABINIT,Quantum
%     We recommend ABINIT,QE and Elk. Default set from installed software.
%   options.dos=1                          Option to compute the vibrational
%     density of states (vDOS) in UserData.DOS
%   options.potentials=string              Basis datasets or pseudopotentials.
%     NWChem: see http://www.nwchem-sw.org/index.php/Release64:AvailableBasisSets
%       e.g. 'ANO-RCC'
%     Dacapo: path to potentials, e.g. /usr/share/dacapo-psp       (DACAPOPATH)
%     Elk:    path to potentials, e.g. /usr/share/elk-lapw/species (ELK_SPECIES_PATH)
%     ABINIT: path to potentials, e.g. /usr/share/abinit/psp/      (ABINIT_PP_PATH)
%       or any of 'NC' 'fhi' 'hgh' 'hgh.sc' 'hgh.k' 'tm' 'paw'
%     QuantumEspresso: path to potentials, e.g. /usr/share/espresso/pseudo
%   options.command='exe'                  Path to calculator executable.
%   options.mpi=scalar                     use multi-cpu, requires e.g. OpenMPI to be installed.
%   options.autoplot=0|1                   when set, automatically evaluates the Model
%     and plot results.
%   options.htmlreport=0|1                 when set, automatically generates a full
%     report on the computation results. Also requests vDOS computation (options.dos=1).
%   options.email=<email>                  when set, sends an email at start and end
%     of computation.
%
% DFT specific options
%   options.kpoints=scalar or [nx ny nz]   Monkhorst-Pack grid, default 3.
%   options.xc=string                      Type of Exchange-Correlation functional to use
%     'LDA','PBE','revPBE','RPBE','PBE0','B3LYP'   for GPAW
%     'PZ', 'PBE','revPBE','RPBE','PW91','VWN'     for Dacapo/Jacapo
%     'LDA','PBE','revPBE','RPBE','GGA'            for ABINIT
%     'LDA','PBE','REVPBE','PBESOL','WC06','AM05'  for ELK
%     'LDA','PBE','RHF','MP2','B3LYP'              for NWChem
%     'LDA','PBE','PW91'                           for VASP
%     Default is 'PBE'.
%   options.raw='name=value, ...'          Additional arguments passed to the ASE
%                                          using the Python syntax.
%   options.nbands=int                     Number of valence bands for which 
%     wavefunctions are being computed. (= # of electrons /2); for a metal, 20% more
%     (minimum 4 more).
%   options.ecut=scalar                    Kinetic energy cutoff for wavefunctions. 
%     Large value improves convergence. Estimate is 200*n_typ_atoms in [eV].
%     Usually one runs at calculations at various ecut to investigate convergence
%   options.diagonalization='dav' or 'cg' or 'rmm-diis'
%     The Davidson method is faster than the conjugate-gradient but uses more
%     memory. The RMM-DIIS method allows parallelization.
%   options.occupations='metal'            for metals ('smearing') help converge
%                       'insulator'        for insulators
%                       'semiconductor'    sets 0 eV FermiDirac smearing
%                       'auto'             for Elk (automatic smearing)
%                       or 0 for semi-conductors
%                       or a value in eV for the distribution width e.g. 0.3
%   options.nsteps=scalar                   Max number of iterations for SCF.
%     Typical: 30. Large value improve convergence.
%   options.toldfe=scalar                  Convergence threshold on the energy for 
%     selfconsistency, e.g. 1e-5 [eV].
%
% Options specific per calculator
%   options.mode='pw','fd', or 'lcao'      GPAW computation mode as Plane-Wave,
%     Finite Difference, or LCAO (linear combination of atomic orbitals). Default is 'fd'.
%   options.iscf='NC','PAW'                Type of SCF cycles (ABINIT) 
%   options.pps = 'fhi' 'hgh' 'hgh.sc' 'hgh.k' 'tm' 'paw' Type of database (ABINIT)
%                 'sv','pv', ... (VASP)
%   options.mixing_beta=scalar             mixing factor for self-consistency
%     default=0.7. use 0.3 to improve convergence (QuantumEspresso)
%   options.mixing_ndim=scalar             number of iterations used in mixing
%     default=8. If you are tight with memory, you may reduce it to 4. A larger
%     value will improve the SCF convergence, and use more memory (QuantumEspresso)
%
% The options can also be entered as strings with 'field=value; ...'.
%
% output (model creation):
%   model: iFunc 4D model S(qh,qk,ql,w) with parameters p.
%
% Once the model has been created, its use requires that axes are given on
% regular qx,qy,qz grids (in rlu along reciprocal axes). The model evaluations
% does not require to recompute the forces, and is very fast. To generate a
% powder 2D S(q,w) you may use: sqw_powder(model)
%     
% Example (model creation and evaluation):
%   s=sqw_phonons('bulk("Cu", "fcc", a=3.6, cubic=True)','EMT','metal','dos');
%   qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);
%   f=iData(s,[],qh,qk,ql,w); scatter3(log(f(1,:, :,:)),'filled');
%   figure; plot(s.UserData.DOS); % plot the DOS, as indicated during model creation
%
%   s=sqw_phonons('bulk("Si", "diamond", a=5.4, cubic=True)','semiconductor');
%
%   s=sqw_phonons([ ifitpath 'Data/POSCAR_Al'],'dos','metal','EMT');
%
% You may look at the following resources to get material structure files:
%   <http://phonondb.mtl.kyoto-u.ac.jp/>
%   <https://www.materialsproject.org/>
%   <http://nomad-repository.eu/cms/>
%
% WARNING: Single intensity and line width parameters are used here.
%   This model is suitable to compute phonon dispersions for e.g solid-
%   state materials.
%   The Atomic Simulation Environment must be installed.
%   The temporary directories (UserData.dir) are not removed.
%   The intensity is currently not computed, only the dispersions.
%
% References: https://en.wikipedia.org/wiki/Phonon
%
% Atomic Simulation Environment
%   S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002
%   <https://wiki.fysik.dtu.dk/ase>. Exists as Debian package 'python-ase'.
%   Repository: http://download.opensuse.org/repositories/home:/dtufys/
% GPAW J. J. Mortensen et al, Phys Rev B, Vol. 71, 035109 (2005). 
%   <http://wiki.fysik.dtu.dk/gpaw>. Exists as Debian package 'gpaw' and 'gpaw-setups' (gpaw-data).
% NWChem M. Valiev et al, Comput. Phys. Commun. 181, 1477 (2010).
%   <http://www.nwchem-sw.org/>. Exists as Debian package 'nwchem' and 'nwchem-data'.
% Elk <http://elk.sourceforge.net>. Exists as Debian package 'elk-lapw'. 
%   The Elk executable should be compiled with e.g. 
%     elk/src/modmain.f90:289 maxsymcrys=1024 or larger
% DACAPO B. Hammer et al, Phys.Rev. B 59, 7413 (1999) 
%   <https://wiki.fysik.dtu.dk/dacapo/>. Exists as Debian package 'dacapo' and 'dacapo-psp'.
% ABINIT  X. Gonze et al, Computer Physics Communications 180, 2582-2615 (2009).
%   <http://www.abinit.org/>. Exists as Debian package 'abinit' and 'abinit-doc' (abinit-data).
%   Install potentials from http://wiki.fysik.dtu.dk/abinit-files/abinit-pseudopotentials-2.tar.gz
%   into e.g. /usr/share/abinit
% PHON D. Alfe, Computer Physics Communications 180,2622-2633 (2009). 
%   <http://chianti.geol.ucl.ac.uk/~dario>.
% Quantum Espresso P. Giannozzi, et al, J.Phys.:Condens.Matter, 21, 395502 (2009).
%   <http://www.quantum-espresso.org/>. 
% VASP G. Kresse and J. Hafner. Phys. Rev. B, 47:558, 1993.
%   <https://www.vasp.at/> Requires a license.
%
% Once the model has been created:
% input:  p: sqw_phonons model parameters (double)
%             p(1)=Amplitude
%             p(2)=Gamma   dispersion DHO half-width in energy [meV]
%             p(3)=Background (constant)
%             p(4)=Temperature of the material [K]
%          or p='guess'
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_cubic_monoatomic, sqw_sine3d, sqw_vaks
%   sqw_powder, <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>

% Units: 1 Ry        = 13.6 eV
%        1 Ha = 2 Ry = 27.2 eV

% best codes: basis https://molmod.ugent.be/deltacodesdft
% code    basis
% QE      SSSP Accuracy     http://materialscloud.org/sssp/
% Elk     PAW+lo
% VASP    PAW 2015 GW
% QE      SSSP Efficiency
% ABINIT  PAW JTH           http://www.abinit.org/downloads/PAW2 req v7.6

persistent status ld_library_path

if ~exist('ld_library_path') || isempty(ld_library_path) || ~ischar(ld_library_path)
  ld_library_path = getenv('LD_LIBRARY_PATH');
end

if ~exist('status') || isempty(status) || ~isstruct(status)
  status = sqw_phonons_requirements;
end

signal = [];

if nargin == 0, configuration = 'gui'; end

options= sqw_phonons_argin(configuration, varargin{:});
if isempty(status.mpirun) && isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
  options.mpi=1;
  disp([ mfilename ': MPI parallelization is not available. Install e.g. OpenMPI first' ]);
end

if isempty(options.calculator)
  % select default calculator according to those installed
  for index={'quantumespresso','abinit','vasp','gpaw','elk','jacapo','nwchem',};
    if ~isempty(status.(index{1})), options.calculator=index{1}; break; end
  end
end

if strcmp(configuration,'gui')
  options.gui=1;
  configuration=[];
end

if nargin == 0 || isempty(configuration)
  configuration = 'bulk("Al", "fcc", a=4.05)';
end

% ==============================================================================
%                               GUI (dialog)
% ==============================================================================
if options.gui
  % pop-up a simple dialog requesting for:
  %  * configuration
  %  * calculator
  %  * metal, insulator, semiconductor
  %  * supercell
  %  * kpoints
  % and will select autoplot, 'dos'
  calcs = 'EMT';
  for index={'gpaw','elk','jacapo','nwchem','abinit','quantumespresso','vasp'};
    if ~isempty(status.(index{1})), calcs = [ calcs ', ' upper(index{1}) ]; end
  end
  NL = sprintf('\n');
  prompt = { [ '{\bf Atom/molecule/system configuration}' NL 'a CIF/PDB/POSCAR/... name or e.g. bulk("Cu", "fcc", a=3.6, cubic=True),' NL  'molecule("H2O"), or nanotube(6, 0, length=4). ' NL 'Documentation at {\color{blue}http://ifit.mccode.org/Models.html}' ], ...
  [ '{\bf Calculator}' NL 'One of ' calcs ], ...
  [ '{\bf Smearing}' NL 'metal, semiconductor, insulator or left empty. Use e.g. 0.3 eV for conductors, or a small value such as 0.01 to help SCF convergence. You may use "auto" with Elk. ' ], ...
  [ '{\bf Cut-off energy for wave-functions}' NL 'Leave as 0 for the default, or specify a cut-off in eV, e.g. 500 eV for fast estimate, 1500 or 2000 eV for accurate results.' ], ...
  [ '{\bf K-Points}' NL 'Monkhorst-Pack grid which determines the K sampling (3-vector). 4 is the minimum for accurate computations, 6 is best. Use 1 or 2 for testing only (faster).' ], ...
   [ '{\bf Supercell}' NL 'The size of the repeated model = system*supercell (3-vector). Should be larger than k-points.' ], ...
   [ '{\bf Other options}' NL 'Such as mpi, nbands, nsteps, xc (default PBE), toldfe, raw' NL 'example: "mpi=4; nsteps=100"' NL 'Documentation at {\color{blue}http://ifit.mccode.org/Models.html}' ] };
  dlg_title = 'iFit: Model: phonons';
  defAns    = {configuration, options.calculator, options.occupations, num2str(options.ecut), ...
    num2str(options.kpoints), num2str(options.supercell), ''};
  num_lines = [ 1 1 1 1 1 1 1 ]';
  op.Resize      = 'on';
  op.WindowStyle = 'normal';   
  op.Interpreter = 'tex';
  answer = inputdlg(prompt, dlg_title, num_lines, defAns, op);
  if isempty(answer), 
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
  options.dos        = 1;
  options.gui        = waitbar(0, [ mfilename ': starting' ], 'Name', [ 'iFit: ' mfilename ' ' configuration ]);
  if strcmpi(options.calculator,'qe') options.calculator='quantumespresso'; end
end % GUI

% ==============================================================================
%                               BUILD MODEL
% ==============================================================================

t=clock();

% handle compatibility fields
if isfield(options,'conv_thr') && ~isfield(options,'toldfe')
    options.toldfe = options.conv_thr; end
    
if strcmp(options.occupations, 'smearing') || strcmp(options.occupations, 'metal') % metals
    options.occupations=0.27; % recommended by SSSP
  elseif strcmp(options.occupations, 'semiconductor')
    options.occupations=0.0001;
  elseif strcmp(options.occupations, 'fixed') || strcmp(options.occupations, 'insulator') % insulators
    options.occupations=-1;
  end

options.configuration = configuration;

% BUILD stage: we call ASE to build the model

pw = pwd; target = options.target;

% handle input configuration: read
if exist(configuration)
  read = sprintf('import ase.io; configuration = "%s"; atoms = ase.io.read(configuration); ', ...
    configuration);
elseif ischar(configuration)
  read = configuration;
  switch strtok(configuration, ' (')
  case 'bulk'
    read = sprintf('from ase.lattice import bulk; atoms = %s; ', configuration);
  case 'molecule'
    read = sprintf('from ase.structure import molecule; atoms = %s; ', configuration);
  case 'nanotube'
    read = sprintf('from ase.structure import nanotube; atoms = %s; ', configuration);
  case 'crystal'
    read = sprintf('from ase.lattice.spacegroup import crystal; atoms = %s; ', configuration);
  end
end

% used for export of structure to usual files+additional information
% the Mat file should be loaded as: atoms=load(fullfile(target, 'atoms.mat'))
read1 = [ ...
    'write("' fullfile(target, 'configuration.png') '", atoms);' ...
    'write("' fullfile(target, 'configuration.eps') '", atoms);' ...
    'write("' fullfile(target, 'configuration.pov') '", atoms);' ...
    'write("' fullfile(target, 'configuration.cif') '", atoms, "cif");' ...
    'write("' fullfile(target, 'configuration.x3d') '", atoms, "x3d");' ...
    'write("' fullfile(target, 'configuration.pdb') '", atoms, "pdb");' ...
    'write("' fullfile(target, 'configuration.html') '", atoms, "html");' ...
    'write("' fullfile(target, 'configuration.etsf') '", atoms, "etsf");' ...
    'write("' fullfile(target, 'configuration_SHELX.res') '", atoms, "res");' ...
    'write("' fullfile(target, 'configuration_VASP') '", atoms, "vasp");' ... 
    'import scipy.io as sio; ' ...
    'sio.savemat("' fullfile(target, 'atoms.mat') '", {' ...
    '"reciprocal_cell": atoms.get_reciprocal_cell(), ' ...
    '"cell": atoms.get_cell(), ' ...
    '"volume": atoms.get_volume(), ' ...
    '"chemical_formula": atoms.get_chemical_formula(), ' ...
    '"masses": atoms.get_masses(), ' ...
    '"positions": atoms.get_positions(), ' ...
    '"atomic_numbers": atoms.get_atomic_numbers() });  ' ...
    ];

% handle the optimizer
if ~isempty(options.optimizer)
  switch lower(options.optimizer)
  case 'lbfgs'  % fast, low on memory
    options.optimizer='LBFGS';
  case 'fire'   % slow
    options.optimizer='FIRE';
  case 'mdmin'  % fast
    options.optimizer='MDMin';
  otherwise
  % case 'bfgs'   % fast
    options.optimizer='BFGS';
  end
  optim = sprintf([ 'from ase.optimize import %s; atoms.set_calculator(calc); ' ...
         'print "Optimising structure with the %s optimizer...\\n"; ' ...
         'dyn = %s(atoms); dyn.run(fmax=0.05); ' ...
         'from ase.io import write; ' ...
         'write("%s", atoms, "cif");' ...
         'write("%s", atoms, "pdb");' ...
         'write("%s", atoms, "vasp");' ...
         ], ...
            options.optimizer, options.optimizer, options.optimizer, ...
            fullfile(options.target, 'optimized.cif'), fullfile(options.target, 'optimized.pdb'), ...
            fullfile(options.target, 'optimized_POSCAR'));
else
  optim = '';
end


% ------------------------------------------------------------------------------
% handle supported calculators: decl and calc
if options.gui && ishandle(options.gui), waitbar(0.05, options.gui, [ mfilename ': configuring ' options.calculator ]); end

if strcmpi(options.calculator, 'QUANTUMESPRESSO') 
  if ~isempty(options.optimizer)
    % specific case for QE. As it is not directly supported by ASE, the miniasation
    % must use an other ASE-compliant code. We try ABINIT, then GPAW.
    if ~isempty(status.abinit)
      [decl, calc] = sqw_phonons_calc(options, status, 'ABINIT');
    elseif ~isempty(status.gpaw)
      [decl, calc] = sqw_phonons_calc(options, status, 'GPAW');
    else
      disp([ mfilename ': WARNING: can not perform optimization with QuantumEspresso.' ])
      disp('    This step requires to use an other DFT code. Install either ABINIT or GPAW.');
      disp('    Skipping.');
    end
    read = [ read optim ];
  end
  % create the POSCAR input file for PHON
  % ASE is installed. We use it to create a proper POSCAR file, then we call sqw_phon (QE)
  poscar = fullfile(options.target,'POSCAR_ASE');
  read = [ read '; from ase.io import write; ' ...
     'write("' poscar '",atoms, "vasp"); ' ...
     read1 ];
end

if isunix, setenv('LD_LIBRARY_PATH',''); end
[decl, calc] = sqw_phonons_calc(options, status, options.calculator, read);
if isunix, setenv('LD_LIBRARY_PATH',ld_library_path); end
  







if ~strcmpi(options.calculator, 'QUANTUMESPRESSO')

  if options.gui && ishandle(options.gui), waitbar(0.10, options.gui, [ mfilename ': writing python script ASE/' options.calculator ]); end

  if strcmpi(options.calculator, 'GPAW')
    % GPAW Bug: gpaw.aseinterface.GPAW does not support pickle export for 'input_parameters'
    sav = sprintf('ph.calc=None\npickle.dump(ph, fid)');
  else
    sav = 'pickle.dump(ph, fid)';
  end
  % start python --------------------------
  script = { ...
    '# python script built by ifit.mccode.org/Models.html sqw_phonons', ...
  [ '# on ' datestr(now) ], ...
    '# E. Farhi, Y. Debab and P. Willendrup, J. Neut. Res., 17 (2013) 5', ...
    '# S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002.', ...
    '#', ...
    '# Computes the dynamical matrix and stores an ase.phonon.Phonons object in a pickle ph.pkl', ...
    '# Launch with: python sqw_phonons_build.py (and wait...)', ...
    decl, ...
    'from ase.phonons import Phonons', ...
    'import pickle', ...
    '# Setup crystal and calculator', ...
    read, ...
    calc, ...
    optim, ...
    '# Phonon calculator', ...
  sprintf('ph = Phonons(atoms, calc, supercell=(%i, %i, %i), delta=0.05)',options.supercell), ...
    'ph.run()', ...
    '# Read forces and assemble the dynamical matrix', ...
    'ph.read(acoustic=True)', ...
    '# save ph', ...
  [ 'fid = open(''' fullfile(target, 'ph.pkl') ''',''wb'')' ], ...
    sav, ...
    'fid.close()' };
  % end   python --------------------------

  % write the script in the target directory
  fid = fopen(fullfile(target,'sqw_phonons_build.py'),'w');
  fprintf(fid, '%s\n', script{:});
  fclose(fid);
  % copy the configuration into the target
  if ~isempty(dir(configuration))
    try
    copyfile(configuration, target);
    end
  end
  
  % call python script with configuration
  if isunix, precmd = 'LD_LIBRARY_PATH= ; '; else precmd=''; end
  read = [ read 'from ase.io import write; ' ...
    read1 ];
  result = '';
  if isunix, setenv('LD_LIBRARY_PATH',''); end
  try
    [st, result] = system([ precmd 'python -c ''' read '''' ]);
  catch
    st = 127;
  end
  if isunix, setenv('LD_LIBRARY_PATH',ld_library_path); end
  if st ~= 0
    disp(read)
    disp(result)
    sqw_phonons_error([ mfilename ': failed read input ' ...
      configuration ], options);
  end
  sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'init', options, calc);
  % get 'atoms' back from python
  if ~isempty(fullfile(target, 'atoms.mat'))
    atoms = load(fullfile(target, 'atoms.mat'));
  else
    atoms = [];
  end

  % call python script with calculator
  cd(target)
  if ~isdeployed && usejava('jvm') && usejava('desktop')
    disp([ mfilename ': creating Phonon/ASE model in <a href="' target '">' target '</a>' ]);
  else
    disp([ mfilename ': creating Phonon/ASE model in ' target ]);
  end
  disp([ '  ' configuration ]);
  disp([ '  ' calc ]);
  if ~isempty(options.optimizer)
    disp([ '  The initial structure will first be optimized using ' options.optimizer ]);
  end

  if options.gui && ishandle(options.gui), waitbar(0.15, options.gui, [ mfilename ': computing DynMatrix ASE/' options.calculator ' (be patient)' ]); end

  result = '';
  if isunix, setenv('LD_LIBRARY_PATH',''); end
  try
    if strcmpi(options.calculator, 'GPAW') && isfield(options,'mpi') ...
      && ~isempty(options.mpi) && options.mpi > 1
      [st, result] = system([ precmd status.mpirun ' -np ' num2str(options.mpi) ' '  status.gpaw ' sqw_phonons_build.py' ]);
    else
      [st, result] = system([ precmd 'python sqw_phonons_build.py' ]);
    end
    disp(result)
  catch
    disp(result)
    if isunix, setenv('LD_LIBRARY_PATH',ld_library_path); end
    sqw_phonons_error([ mfilename ': failed calling ASE with script ' ...
      fullfile(target,'sqw_phonons_build.py') ], options);
  end
  cd(pw)

  % then read the pickle file to store it into the model
  if options.gui && ishandle(options.gui), waitbar(0.7, options.gui, [ mfilename ': building model' ]); end
  try
    signal.UserData.ph_ase = fileread(fullfile(target, 'ph.pkl')); % binary
  catch
    setenv('LD_LIBRARY_PATH',ld_library_path);
    sqw_phonons_error([ mfilename ': ' options.calculator ' failed. Temporary files and Log are in ' target ], options)
  end
  if ~isempty(dir(configuration))
    [dummy, signal.UserData.input]= fileparts(configuration);
  else
    signal.UserData.input = configuration;
  end

  signal.Name           = [ 'Sqw_ASE_' signal.UserData.input ' Phonon/ASE/' options.calculator ' DHO [' mfilename ']' ];

  signal.Description    = [ 'S(q,w) 3D dispersion Phonon/ASE/' options.calculator ' with DHO line shape. ' configuration ];

  signal.Parameters     = {  ...
    'Amplitude' ...
    'Gamma Damped Harmonic Oscillator width in energy [meV]' ...
    'Background' ...
    'Temperature [K]' ...
     };
    
  signal.Dimension      = 4;         % dimensionality of input space (axes) and result

  signal.Guess = [ 1 .1 0 10 ];

  if ~isempty(dir(configuration))
    signal.UserData.configuration = fileread(configuration);
  else
    signal.UserData.configuration = configuration;
  end
  if isfield(options,'dos') && options.dos, signal.UserData.DOS=[]; end
  signal.UserData.dir           = target;
  signal.UserData.options       = options;
  signal.UserData.calc          = calc;
  signal.UserData.atoms         = atoms;

  % EVAL stage: we call ASE to build the model. ASE does not support single HKL location. 
  % For this case we duplicate xyz and then get only the 1st FREQ line

  signal.Expression = { ...
    '% check if directory and phonon pickle is here', ...
  [ 'pw = pwd; target = this.UserData.dir;' ], ...
    'if ~isdir(target), target = tempname; mkdir(target); this.UserData.dir=target; end', ...
  [ 'if isempty(dir(fullfile(target, ''ph.pkl'')))' ], ...
    '  fid=fopen(fullfile(target, ''ph.pkl''), ''w'');', ...
    '  if fid==-1, error([ ''model '' this.Name '' '' this.Tag '' could not write ph.pkl into '' target ]); end', ...
    '  fprintf(fid, ''%s\n'', this.UserData.ph_ase);', ...
    '  fclose(fid);', ...
    'end', ...
  [ '  fid=fopen(fullfile(target,''sqw_phonons_eval.py''),''w'');' ], ...
  [ '  fprintf(fid, ''# Script file for Python/ASE to compute the modes from ' configuration ' in %s\n'', target);' ], ...
    '  fprintf(fid, ''#   ASE: S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002\n'');', ...
    '  fprintf(fid, ''#   <https://wiki.fysik.dtu.dk/ase>. Built by ifit.mccode.org/Models.html sqw_phonons'');', ...
    '  fprintf(fid, ''from ase.phonons import Phonons\n'');', ...
    '  fprintf(fid, ''import numpy\n'');', ...
    '  fprintf(fid, ''import pickle\n'');', ...
    '  fprintf(fid, ''# restore Phonon model\n'');', ...
    '  fprintf(fid, ''fid = open(''''ph.pkl'''', ''''rb'''')\n'');', ...
    '  fprintf(fid, ''ph = pickle.load(fid)\n'');', ...
    '  fprintf(fid, ''fid.close()\n'');', ...
    '  fprintf(fid, ''# read HKL locations\n'');', ...
    '  fprintf(fid, ''HKL = numpy.loadtxt(''''HKL.txt'''')\n'');', ...
    '  fprintf(fid, ''# compute the spectrum\n'');', ...
    '  fprintf(fid, ''omega_kn = 1000 * ph.band_structure(HKL)\n'');', ...
    '  fprintf(fid, ''# save the result in FREQ\n'');', ...
    '  fprintf(fid, ''numpy.savetxt(''''FREQ'''', omega_kn)\n'');', ...
    'if isfield(this.UserData, ''DOS'') && isempty(this.UserData.DOS)', ...
    '  fprintf(fid, ''# Calculate phonon DOS\n'');', ...
    '  fprintf(fid, ''omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=5e-4)\n'');', ...
    '  fprintf(fid, ''omega_e *= 1000\n'');', ...
    '  fprintf(fid, ''# save the result in DOS\n'');', ...
    '  fprintf(fid, ''numpy.savetxt(''''DOS_w'''',omega_e)\n'');', ...
    '  fprintf(fid, ''numpy.savetxt(''''DOS'''',  dos_e)\n'');', ...
    'end', ...
    '  fprintf(fid, ''exit()\n'');', ...
    '  fclose(fid);', ...
    '  sz0 = size(t);', ...
    '  if ndims(x) == 4, x=squeeze(x(:,:,:,1)); y=squeeze(y(:,:,:,1)); z=squeeze(z(:,:,:,1)); t=squeeze(t(1,1,1,:)); end',...
    'try', ...
    '  cd(target);', ...
    '  if all(cellfun(@isscalar,{x y z t})), HKL = [ x y z ; x y z ];', ...
    '  else HKL = [ x(:) y(:) z(:) ]; end', ...
    '  save -ascii HKL.txt HKL', ...
  [ '  [status,result] = system(''' precmd 'python sqw_phonons_eval.py'');' ], ...
    '  clear HKL', ...
    '  delete(''HKL.txt'')', ...
    '  % import FREQ', ...
    '  FREQ=load(''FREQ'',''-ascii''); % in meV', ...
    '  delete(''FREQ'')', ...
    'catch ME; disp([ ''model '' this.Name '' '' this.Tag '' could not run Python/ASE from '' target ]); disp(getReport(ME))', ...
    'end', ...
    'try', ...
    'if isfield(this.UserData, ''DOS'') && isempty(this.UserData.DOS)', ...
    '  DOS = load(''DOS'',''-ascii''); DOS_w = load(''DOS_w'',''-ascii''); DOS=iData(DOS_w,DOS./sum(DOS));', ...
    '  DOS.Title = [ ''DOS '' strtok(this.Name) ]; xlabel(DOS,''DOS Energy [meV]'');', ...
    '  DOS.Error=0; this.UserData.DOS=DOS;', ...
    'end', ...
    'catch ME; disp([ ''model '' this.Name '' '' this.Tag '' could not get vDOS Python/ASE from '' target ]); disp(getReport(ME));', ...
    'end', ...
    '  cd(pw);', ...
    '  % multiply all frequencies(columns, meV) by a DHO/meV', ...
    '  Amplitude = p(1); Gamma=p(2); Bkg = p(3); T=p(4);', ...
    '% apply DHO on each',...
    '  if T<=0, T=300; end', ...
    '  if all(cellfun(@isscalar,{x y z t})), FREQ=FREQ(1,:); end', ...
    '  if all(cellfun(@isvector,{x y z t}))', ...
    '  w=t(:); else w=t(:) * ones(1,size(FREQ,1)); end', ...
    '  signal=zeros(size(w));', ...
    'for index=1:size(FREQ,2)', ...
    '% transform w and w0 to same size', ...
    '  if all(cellfun(@isvector,{x y z t}))', ...
    '  w0 =FREQ(:,index); else w0= ones(numel(t),1) * FREQ(:,index)''; end', ...
    '  toadd = Amplitude*Gamma *w0.^2.* (1+1./(exp(abs(w)/T)-1)) ./ ((w.^2-w0.^2).^2+(Gamma*w).^2);', ...
    '  signal = signal +toadd;', ...
    'end', ...
    'signal = reshape(signal'',sz0);'};

  signal = iFunc(signal);
end % other calculators than QE

signal.UserData.duration = etime(clock, t);
options.duration = signal.UserData.duration;
setenv('LD_LIBRARY_PATH',ld_library_path);

% when model is successfully built, display citations
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '<a href="matlab:doc(''' mfilename ''')">' mfilename '</a>: Model ' configuration ' built using ' options.calculator ])
else
  disp([ mfilename ': Model ' configuration ' built using ' options.calculator ])
end
if isdeployed || ~usejava('jvm') || ~usejava('desktop')
  disp([ '  in ' options.target ]);
else
  disp([ '  in <a href="' options.target '">' options.target '</a>' ]);
end


if isfield(options, 'dos') && options.dos && ~strcmpi(options.calculator, 'QUANTUMESPRESSO')
  disp('INFO: The vibrational density of states (vDOS) will be computed at first model evaluation.');
end
disp([ 'Time elapsed=' num2str(signal.UserData.duration) ' [s]. Please cite:' ])
cite = { ' * Atomic Simulation Environment', ...
  '           S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002.', ...
  ' * iFit:   E. Farhi et al, J. Neut. Res., 17 (2013) 5.' };
switch upper(options.calculator)
case 'GPAW'
cite{end+1} = ' * GPAW:   J. J. Mortensen et al, Phys Rev B, Vol. 71, 035109 (2005).';
case 'NWCHEM'
cite{end+1} = ' * NWChem: M. Valiev et al, Comput. Phys. Commun. 181, 1477 (2010).';
case 'ELK'
cite{end+1} = ' * ELK:    http://elk.sourceforge.net';
case {'DACAPO','JACAPO'}
cite{end+1} = ' * DACAPO: B. Hammer et al, Phys. Rev. B 59, 7413 (1999).';
case 'ABINIT'
cite{end+1} = ' * ABINIT: X. Gonze et al, Computer Physics Communications 180, 2582-2615 (2009).';
case 'EMT'
cite{end+1} = ' * EMT:    K.W. Jacobsen et al, Surf. Sci. 366, 394-402 (1996).';
case 'QUANTUMESPRESSO'
cite{end+1} = ' * PHON:   D. Alfe, Computer Physics Communications 180,2622-2633 (2009).';
cite{end+1} = ' * Quantum Espresso: P. Giannozzi, et al J.Phys.:Condens.Matter, 21, 395502 (2009).';
case 'VASP'
cite{end+1} = ' * VASP:   G. Kresse and J. Hafner. Phys. Rev. B, 47:558, 1993.';
end
fprintf(1, '%s\n', cite{:});
disp('You can now evaluate the model using e.g.:')
disp('    qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);');
disp('    f=iData(s,[],qh,qk,ql,w); plot3(log(f(1,:, :,:)));');
disp(' ');

% save the Model as a Matlab object
Phonons_Model = signal;
builtin('save', fullfile(options.target, 'Phonons_Model.mat'), 'Phonons_Model');

sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'done', options, cite);

% handle autoplot option
if options.autoplot
  if options.gui && ishandle(options.gui), waitbar(0.75, options.gui, [ mfilename ': plotting phonons and DOS' ]); end
  [f, signal] = sqw_phonons_plot(signal);
  if options.gui && ishandle(options.gui), delete(options.gui); end
else
  f = [];
end

sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'plot', options, f, signal);
sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'end', options, f, signal);

% ------------------------------------------------------------------------------
function [f, signal] = sqw_phonons_plot(signal)
  if isempty(signal) || ~isfield(signal.UserData, 'options') , return; end
  
  options = signal.UserData.options;
  disp([ mfilename ': Model ' options.configuration ' plotting phonons.' ])
  qh=linspace(0.01,1.5,70);qk=qh; ql=qh; w=linspace(0.01,150,151);
  f=iData(signal,[],qh,qk,ql,w);
  g=log(f(1,:, :,:)); 
  
  fig1=figure; 
  plot3(g); axis tight;
  view([38 26]);
  if isfield(options, 'dos') && options.dos && ~isempty(signal.UserData.DOS)
    fig2=figure;
    try
    plot(signal.UserData.DOS); % plot the DOS, as indicated during model creation
    catch
    close(fig2)
    end
  end
  drawnow

% ------------------------------------------------------------------------------
function sqw_phonons_error(message, options)

if options.gui && ishandle(options.gui)
  delete(options.gui);
  errordlg(message, [ 'iFit: ' mfilename ' ' options.configuration ' FAILED' ]);
end
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '<a href="matlab:doc(''' mfilename ''')">help ' mfilename '</a> (click here to get help)' ])
end
sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'error', options, message);
error(message);

% ------------------------------------------------------------------------------



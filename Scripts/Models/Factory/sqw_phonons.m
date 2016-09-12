function signal=sqw_phonons(configuration, varargin)
% model=sqw_phonons(configuration, calculator, ..., options)
%
%   iFunc/sqw_phonons: computes phonon dispersions using the ASE.
%   A model which computes phonon dispersions from the Hellmann-Feynman forces acting 
%     between atoms. The input argument is any configuration file describing the
%     material, e.g. CIF, PDB, POSCAR, ... supported by ASE.
%   This models can compute the coherent inelastic phonons dispersions for
%     any crystalline (powder or single crystal) material, in the harmonic 
%     approxiimation.
%   The phonon spectra is computed using one of the calculator supported by the
%   Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase>.
%
%   Supported calculators are:
%     ABINIT    Plane-wave pseudopotential code
%     Dacapo    Plane-wave ultra-soft pseudopotential code
%     ELK       Full Potential LAPW code
%     EMT       Effective Medium Theory calculator (Al,Cu,Ag,Au,Ni,Pd,Pt,H,C,N,O)
%     GPAW      Real-space/plane-wave/LCAO PAW code
%     NWChem    Gaussian based electronic structure code
%     QuantumEspresso_ASE Plane-wave pseudopotential code (with ASE/QE-util)
%     QuantumEspresso     Plane-wave pseudopotential code (with PHON)
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
%   * QuantumEspresso/PHON is the fastest, with excellent parallelization and accuracy.
%   * QuantumEspresso/ASE is excellent.
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
%   options.optimizer                      when set, performs a geometry optimization
%     before computing the forces. This option can be set to 
%       'BFGS' or 'LBFGS' (low memory BFGS)
%       'MDMin' or 'FIRE'
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
%     ABINIT paw/pawxml potentials require larger ecut (1000 eV) and nsteps (100).
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
%     Typical: 30. Large value improves convergence.
%   options.toldfe=scalar                  Convergence threshold on the energy for 
%     selfconsistency, e.g. 1e-5 [eV].
%
% Options specific per calculator
%   options.mode='pw','fd', or 'lcao'      GPAW computation mode as Plane-Wave,
%     Finite Difference, or LCAO (linear combination of atomic orbitals). Default is 'pw'.
%   options.iscf='NC','PAW'                Type of SCF cycles (ABINIT) 
%   options.pps = 'fhi' 'hgh' 'hgh.sc' 'hgh.k' 'tm' 'paw' 'pawxml' Type of database (ABINIT)
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
% (c) E.Farhi, ILL. License: EUPL.

% Units: 1 Ry        = 13.6 eV
%        1 Ha = 2 Ry = 27.2 eV
%        1 meV = 241.8 GHz = 11.604 K = 0.0965 kJ/mol = 8.0657 cm-1 

% best codes: basis https://molmod.ugent.be/deltacodesdft
% code    basis
% QE      SSSP Accuracy     http://materialscloud.org/sssp/
% Elk     PAW+lo
% VASP    PAW 2015 GW
% QE      SSSP Efficiency
% ABINIT  PAW JTH           http://www.abinit.org/downloads/PAW2 req v7.6
%                           https://www.physics.rutgers.edu/gbrv/

% python scripts are writen for ASE:
%   sqw_phonons_check.py:     read the initial structure and write standard format
%   sqw_phonons_optimize.py:  optimize the structure and return a new configuration
%   sqw_phonons_build.py:     compute the FORCES
%   sqw_phonons_eval.py:      evaluate the Model at given HKLE locations

persistent status ld_library_path

if ~exist('ld_library_path') || isempty(ld_library_path) || ~ischar(ld_library_path)
  ld_library_path = getenv('LD_LIBRARY_PATH');
end

if ~exist('status') || isempty(status) || ~isstruct(status)
  status = sqw_phonons_requirements;
end

signal = [];

if nargin == 0, configuration = 'gui'; varargin{1} = 'emt'; 
elseif strcmp(configuration, 'identify')
  signal = sqw_phonons('defaults');
  signal.Name = [ 'Sqw_phonon DHO [' mfilename ']' ];
  return;
end

options= sqw_phonons_argin(configuration, varargin{:});
if isempty(status.mpirun) && isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
  options.mpi=1;
  disp([ mfilename ': MPI parallelization is not available. Install e.g. OpenMPI first' ]);
end

if isempty(options.calculator)
  % select default calculator according to those installed
  for index={'quantumespresso','quantumespresso_ase',  ...
      'abinit','vasp','gpaw','elk','jacapo','nwchem'};
    if ~isempty(status.(index{1})), options.calculator=index{1}; break; end
  end
end

if strcmp(configuration,'gui')
  options.gui='init';
  configuration=[];
end

if strcmp(configuration, 'defaults')
  options.calculator = 'emt';
end
if nargin == 0 || isempty(configuration) || strcmp(configuration, 'defaults')
  configuration = 'bulk("Al", "fcc", a=4.05)';
end


% ==============================================================================
%                               GUI (dialog)
% ==============================================================================

if ~isempty(options.gui) && ~any(isnan(options.gui))
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
  options.gui        = true;
  if strcmpi(options.calculator,'qe_ase') options.calculator='quantumespresso_ase'; 
  elseif strcmpi(options.calculator,'qe')         options.calculator='quantumespresso'; 
  end
end % GUI

% make sure we find the 'configuration' also from ifitpath
if isempty(dir(configuration)) && ~isempty(dir(fullfile(ifitpath,'Data',configuration)))
  configuration = fullfile(ifitpath,'Data',configuration);
end

% ==============================================================================
%                               BUILD MODEL (read initial structure)
% ==============================================================================

t = clock();

if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
else           precmd = ''; end

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
  read = sprintf('import ase.io\nconfiguration = "%s"\natoms = ase.io.read(configuration)\n', ...
    configuration);
elseif ischar(configuration)
  read = configuration;
  % ASE has changed some of the modules hierarchy from 3.9 to 3.10+
  switch strtok(configuration, ' (')
  case 'bulk'
    % ASE 3.9 and 3.10+
    read = sprintf('from ase.lattice import bulk\natoms = %s\n', configuration);
  case 'molecule'
    if status.ase <= 3.9
      % ASE 3.9
      read = sprintf('from ase.structure import molecule\natoms = %s\n', configuration);
    else
      % ASE 3.10+: from ase.build import molecule
      read = sprintf('from ase.build import molecule\natoms = %s\n', configuration);
    end
  case 'nanotube'
    if status.ase <= 3.9
      % ASE 3.9
      read = sprintf('from ase.structure import nanotube\natoms = %s\n', configuration);
    else
      % ASE 3.10+: from ase.build import nanotube
      read = sprintf('from ase.build import nanotube\natoms = %s\n', configuration);
    end
  case 'crystal'
    if status.ase <= 3.9
      % ASE 3.9
      read = sprintf('from ase.lattice.spacegroup import crystal\natoms = %s\n', configuration);
    else
      % ASE 3.10+: from ase.spacegroup import crystal
      read = sprintf('from ase.spacegroup import crystal\natoms = %s\n', configuration);
    end
  end
end

options.script_create_atoms = sprintf([ ...
  '# python script built by ifit.mccode.org/Models.html sqw_phonons\n' ...
  '# on ' datestr(now) '\n' ...
  '# E. Farhi, Y. Debab and P. Willendrup, J. Neut. Res., 17 (2013) 5\n' ...
  '# S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002.\n' ...
  '#\n' ...
  '# read initial material structure and save atoms as a pickle\n\n' ...
  read ...
  'import pickle\n' ...
  'fid = open("' fullfile(target, 'atoms.pkl') '","wb")\n' ...
  'pickle.dump(atoms, fid)\n' ...
  'fid.close()' ...
  '# export atoms into usual formats\n' ...
  'print "Exporting structure to usual formats...\\n"\n' ...
  'from ase.io import write\n' ...
  'write("' fullfile(target, 'configuration.png') '", atoms)\n' ...
  'write("' fullfile(target, 'configuration.eps') '", atoms)\n' ...
  'write("' fullfile(target, 'configuration.pov') '", atoms)\n' ...
  'write("' fullfile(target, 'configuration.cif') '", atoms, "cif")\n' ...
  'write("' fullfile(target, 'configuration.x3d') '", atoms, "x3d")\n' ...
  'write("' fullfile(target, 'configuration.pdb') '", atoms, "pdb")\n' ...
  'write("' fullfile(target, 'configuration.html') '", atoms, "html")\n' ...
  'write("' fullfile(target, 'configuration.etsf') '", atoms, "etsf")\n' ...
  'write("' fullfile(target, 'configuration_SHELX.res') '", atoms, "res")\n' ...
  'write("' fullfile(target, 'configuration_VASP') '", atoms, "vasp")\n' ... 
  'import scipy.io as sio\n' ...
  'print "Exporting structure properties...\\n"\n' ...
  'properties = {\n' ...
  '  "reciprocal_cell" : atoms.get_reciprocal_cell(), \n' ...
  '  "cell"            : atoms.get_cell(), \n' ...
  '  "volume"          : atoms.get_volume(), \n' ...
  '  "chemical_formula": atoms.get_chemical_formula(), \n' ...
  '  "chemical_symbols": atoms.get_chemical_symbols(), \n' ...
  '  "cell_scaled_unit": atoms.get_scaled_positions(), \n' ...
  '  "masses"          : atoms.get_masses(), \n' ...
  '  "positions"       : atoms.get_positions(), \n' ...
  '  "atomic_numbers"  : atoms.get_atomic_numbers() }\n' ...
  'try:\n' ...
  '  from pyspglib import spglib\n' ...
  '  properties["spacegroup"] = spglib.get_spacegroup(atoms)\n' ...
  'except:\n' ...
  '  pass\n' ...
  '# export properties as pickle\n' ...
  'fid = open("' fullfile(target, 'properties.pkl') '","wb")\n' ...
  'pickle.dump(properties, fid)\n' ...
  'fid.close()\n' ...
  '# export properties as MAT\n' ...
  'sio.savemat("' fullfile(target, 'properties.mat') '", properties)\n' ...
  ]);

% we create a python script to read/check/export the initial structure
fid = fopen(fullfile(target,'sqw_phonons_check.py'),'w');
if fid == -1
  sqw_phonons_error([ mfilename ': failed create material structure reader for ' ...
    configuration ' (sqw_phonons_check.py)' ], options);
  return
end
fprintf(fid, '%s\n', options.script_create_atoms);
fclose(fid);

% display message at start
disp(' ')
disp([ mfilename ': starting phonons computation [' datestr(now) ']' ])
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '  directory  = <a href="' target '">' target '</a>' ]);
else
  disp([ '  directory  = ' target ]);
end
disp(sprintf([ '  material   = ' configuration ]));
disp(sprintf([ '  calculator = ' options.calculator ]));

% ------------------------------------------------------------------------------
% call python: read/export initial configuration
result = '';
try
  [st, result] = system([ precmd 'python ' fullfile(target,'sqw_phonons_check.py') ]);
end
% we test if the pickle file could be written. This way even if the export/save 
% properties fail, we can proceed.
disp(result)
if isempty(dir(fullfile(target, 'atoms.pkl')))  % FATAL
  sqw_phonons_error([ mfilename ': failed read input ' configuration ], options);
  return
end

% sqw_phonons_htmlreport(fullfile(options.target, 'sqw_phonons.html'), 'init', options, calc);

% ==============================================================================
%                               BUILD MODEL (get calculator)
% ==============================================================================

% get the calculator
% for QE, this triggers computation with sqw_phon()
[decl, calc, signal] = sqw_phonons_calc(options, status, options.calculator, read);
if isempty(decl) && isempty(signal), return; end

% ==============================================================================
%                               BUILD MODEL (optimize or pass)
% ==============================================================================

% handle the optimizer
% requires 'atoms' and 'calc'. provides new 'atoms'.

% init calculator
if strcmpi(options.calculator, 'GPAW')
  % GPAW Bug: gpaw.aseinterface.GPAW does not support pickle export for 'input_parameters'
  sav = sprintf('  atoms.set_calculator(None)\n');
else
  sav = '';
end

if ~isempty(options.optimizer) && ...
  (~strcmpi(options.calculator, 'QUANTUMESPRESSO') || ...
    strcmpi(options.calculator, 'QUANTUMESPRESSO_ASE'))
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
  options.script_optimize = sprintf([ ...
    '# python script built by ifit.mccode.org/Models.html sqw_phonons\n' ...
    '# on ' datestr(now) '\n' ...
    '# E. Farhi, Y. Debab and P. Willendrup, J. Neut. Res., 17 (2013) 5\n' ...
    '# S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002.\n' ...
    '#\n' ...
    '# optimize material structure and update atoms as a pickle\n\n' ...
    'import pickle\n' ...
    'from ase.optimize import ' options.optimizer '\n' ...
    'fid = open("' fullfile(target, 'atoms.pkl') '","rb")\n' ...
    'atoms = pickle.load(fid)\n' ...
    'fid.close()\n' ...
    decl '\n' ...
    calc '\n' ...
    'print "Optimising structure with the ' options.optimizer ' optimizer...\\n"\n' ...
    'try:\n' ...
    '  atoms.set_calculator(calc)\n' ...
    '  dyn = ' options.optimizer '(atoms)\n' ...
    '  dyn.run(fmax=0.05)\n' ...
    '# update pickle and export (overwrite initial structure)\n' ...
    sav ...
    '  fid = open("' fullfile(target, 'atoms.pkl') '","wb")\n' ...
    '  pickle.dump(atoms, fid)\n' ...
    '  fid.close()\n' ...
    '  from ase.io import write\n' ...
    '  write("' fullfile(target, 'optimized.cif') '", atoms, "cif")\n' ...
    '  write("' fullfile(target, 'optimized.pdb') '", atoms, "pdb")\n' ...
    '  write("' fullfile(target, 'optimized_POSCAR') '", atoms, "vasp")\n' ...
    'except:\n' ...
    '  print "Optimisation failed. Proceeding with initial structure\\n"\n' ...
    ]);
    
  % we create a python script to optimize the initial structure
  fid = fopen(fullfile(target,'sqw_phonons_optimize.py'),'w');
  if fid ~= -1
    fprintf(fid, '%s\n', options.script_optimize);
    fclose(fid);
    % call python: optimize initial configuration
    % ------------------------------------------------------------------------------
    result = '';
    disp([ mfilename ': optimizing material structure.' ]);
    try
      [st, result] = system([ precmd 'python ' fullfile(target,'sqw_phonons_optimize.py') ]);
    catch
      st = 127;
    end
    disp(result)
    % was there an error ? if so, continue (we still use the initial atoms.pkl)
    if isempty(dir(fullfile(target, 'atoms.pkl'))) || st ~= 0
      disp([ mfilename ': WARNING: failed optimize material structure in ' target ' (sqw_phonons_optimize.py). Ignoring.' ]);
    else
      % sqw_phonons_htmlreport(fullfile(options.target, 'sqw_phonons.html'), 'optimize', options, calc);
    end

  end
else
  options.script_optimize = '';
end



% ==============================================================================
%                               COMPUTE MODEL (force computation)
% ==============================================================================

% the QE case has been done prior to optimize when calling sqw_phonons_calc()

if ~strcmpi(options.calculator, 'QUANTUMESPRESSO') || strcmpi(options.calculator, 'QUANTUMESPRESSO_ASE')

  % init calculator
  if strcmpi(options.calculator, 'GPAW')
    % GPAW Bug: gpaw.aseinterface.GPAW does not support pickle export for 'input_parameters'
    sav = sprintf('ph.calc=None\natoms.calc=None\nph.atoms.calc=None\n');
  else
    sav = '';
  end

  % start python --------------------------
  options.script_get_forces = [ ...
    '# python script built by ifit.mccode.org/Models.html sqw_phonons\n', ...
    '# on ' datestr(now) '\n' ...
    '# E. Farhi, Y. Debab and P. Willendrup, J. Neut. Res., 17 (2013) 5\n', ...
    '# S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002.\n', ...
    '#\n', ...
    '# Computes the Hellmann-Feynman forces and stores an ase.phonon.Phonons object in a pickle\n', ...
    '# Launch with: python sqw_phonons_build.py (and wait...)\n', ...
    'from ase.phonons import Phonons\n', ...
    'import pickle\n', ...
    'import numpy\n', ...
    'import scipy.io as sio\n', ...
    'from os import chdir\n' ...
    'chdir("' target '")\n' ...
    '# Get the crystal and calculator\n', ...
    'fid = open("' fullfile(target, 'atoms.pkl') '","rb")\n' ...
    'atoms = pickle.load(fid)\n' ...
    'fid.close()\n' ...
    decl '\n', ...
    calc '\n', ...
    'atoms.set_calculator(calc)\n' ...
    '# Phonon calculator\n', ...
    sprintf('ph = Phonons(atoms, calc, supercell=(%i, %i, %i), delta=0.05)\n',options.supercell), ...
    'print "Computing Forces\\n"\n', ...
    'ph.run()\n', ...
    '# Read forces and assemble the dynamical matrix\n', ...
    'ph.read(acoustic=True, cutoff=None) # cutoff in Angs\n', ...
    '# save FORCES and phonon object as a pickle\n', ...
    'sio.savemat("' fullfile(target, 'FORCES.mat') '", { "FORCES": ph.C_N , "delta": ph.delta })\n', ...
    'fid = open("' fullfile(target, 'phonon.pkl') '","wb")\n' , ...
    'calc = ph.calc\n', ...
    sav, ...
    'pickle.dump(ph, fid)\n', ...
    'fid.close()\n', ...
    '# additional information\n', ...
    'atoms.set_calculator(calc) # reset calculator as we may have cleaned it for the pickle\n', ...
    'print "Computing properties\\n"\n', ...
    'from ase.vibrations import Vibrations\n', ...
    'try:    magnetic_moment    = atoms.get_magnetic_moment()\n', ...
    'except: magnetic_moment    = None\n', ...
    'try:    kinetic_energy     = atoms.get_kinetic_energy()\n', ... 
    'except: kinetic_energy     = None\n', ...
    'try:    potential_energy   = atoms.get_potential_energy()\n',... 
    'except: potential_energy   = None\n', ...
    'try:    stress             = atoms.get_stress()\n', ... 
    'except: stress             = None\n', ...
    'try:    total_energy       = atoms.get_total_energy()\n', ...
    'except: total_energy       = None\n', ...
    'try:    angular_momentum   = atoms.get_angular_momentum()\n', ... '
    'except: angular_momentum   = None\n', ...
    'try:    charges            = atoms.get_charges()\n', ...
    'except: charges            = None\n', ...
    'try:    dipole_moment      = atoms.get_dipole_moment()\n', ... 
    'except: dipole_moment      = None\n', ...
    'try:    momenta            = atoms.get_momenta()\n', ... 
    'except: momenta            = None\n', ...
    'try:    moments_of_inertia = atoms.get_moments_of_inertia()\n', ...
    'except: moments_of_inertia = None\n', ...
    'try:    center_of_mass     = atoms.get_center_of_mass()\n', ...
    'except: center_of_mass     = None\n', ...
    '# get the previous properties from the init phase\n' ...
    'try:\n' ...
    '  fid = open("' fullfile(target, 'properties.pkl') '","rb")\n' ...
    '  properties = pickle.load(fid)\n' ...
    '  fid.close()\n' ...
    'except:\n' ...
    '  properties = dict()\n' ...
    'try:\n' ...
    '  vib = Vibrations(atoms)\n', ...
    '  print "Computing molecular vibrations\\n"\n', ...
    '  vib.run()\n', ...
    '  vib.summary()\n' ...
    '  properties["zero_point_energy"] = vib.get_zero_point_energy()\n' ...
    '  properties["vibrational_energies"]=vib.get_energies()\n' ...
    'except:\n' ...
    '  pass\n' ...
    'properties["magnetic_moment"]  = magnetic_moment\n' ...
    'properties["kinetic_energy"]   = kinetic_energy\n' ...
    'properties["potential_energy"] = potential_energy\n' ...
    'properties["stress"]           = stress\n' ...
    'properties["momenta"]          = momenta\n' ...
    'properties["total_energy"]     = total_energy\n' ...
    'properties["angular_momentum"] = angular_momentum\n' ...
    'properties["charges"]          = charges\n' ...
    'properties["dipole_moment"]    = dipole_moment\n' ...
    'properties["moments_of_inertia"]= moments_of_inertia\n' ...
    'properties["center_of_mass"]   = center_of_mass\n' ...
    '# remove None values\n' ...
    'properties = {k: v for k, v in properties.items() if v is not None}\n' ...
    '# export properties as pickle\n' ...
    'fid = open("' fullfile(target, 'properties.pkl') '","wb")\n' ...
    'pickle.dump(properties, fid)\n' ...
    'fid.close()\n' ...
    '# export properties as MAT\n' ...
    'sio.savemat("' fullfile(target, 'properties.mat') '", properties)\n' ...
  ];
  % end   python --------------------------

  % write the script in the target directory
  fid = fopen(fullfile(target,'sqw_phonons_forces.py'),'w');
  fprintf(fid, options.script_get_forces);
  fclose(fid);
  % copy the configuration into the target
  if ~isempty(dir(configuration))
    try
    copyfile(configuration, target);
    end
  end
  
  % call python script with calculator
  disp([ mfilename ': computing Hellmann-Feynman forces and creating Phonon/ASE model.' ]);

  result = '';
  try
    if strcmpi(options.calculator, 'GPAW') && isfield(options,'mpi') ...
      && ~isempty(options.mpi) && options.mpi > 1
      [st, result] = system([ precmd status.mpirun ' -np ' num2str(options.mpi) ' '  status.gpaw ' ' fullfile(target,'sqw_phonons_forces.py') ]);
    else
      [st, result] = system([ precmd 'python ' fullfile(target,'sqw_phonons_forces.py') ]);
    end
    disp(result)
  catch
    disp(result)
    sqw_phonons_error([ mfilename ': failed calling ASE with script ' ...
      fullfile(target,'sqw_phonons_build.py') ], options);
    return
  end

  % then read the pickle file to store it into the model
  try
    signal.UserData.phonons = fileread(fullfile(target, 'phonon.pkl')); % binary
  catch
    sqw_phonons_error([ mfilename ': ' options.calculator ' failed. Temporary files and Log are in ' target ], options)
    return
  end
  if ~isempty(dir(configuration))
    [dummy, signal.UserData.input]= fileparts(configuration);
  else
    signal.UserData.input = configuration;
  end
  
  signal.Name           = [ 'Sqw_Phonon_' signal.UserData.input ' ' options.calculator ' DHO [' mfilename ']' ];

  signal.Description    = [ 'S(q,w) dispersion(HKL) Phonon/ASE/' options.calculator ' with DHO line shape. ' configuration ];

  signal.Parameters     = {  ...
    'Amplitude' ...
    'Gamma Damped Harmonic Oscillator width in energy [meV]' ...
    'Background' ...
    'Temperature [K]' ...
     };
    
  signal.Dimension      = 4;         % dimensionality of input space (axes) and result

  signal.Guess = [ 1 .1 0 10 ];
  
  if ~isempty(fullfile(target, 'FORCES.mat'))
    signal.UserData.FORCES = load(fullfile(target, 'FORCES.mat'));
  end
  
  % get 'atoms' back from python
  if ~isempty(fullfile(target, 'properties.mat'))
    properties = load(fullfile(target, 'properties.mat'));
  else
    properties = [];
  end
  
  if ~isempty(dir(configuration))
    signal.UserData.configuration = fileread(configuration);
  else
    signal.UserData.configuration = configuration;
  end
  if isfield(options,'dos') && options.dos, signal.UserData.DOS=[]; end
  signal.UserData.dir           = target;
  signal.UserData.options       = options;
  signal.UserData.calc          = calc;
  signal.UserData.properties    = properties;

  % EVAL stage: we call ASE to build the model. ASE does not support single HKL location. 
  % For this case we duplicate xyz and then get only the 1st FREQ line

  signal.Expression = { ...
    '% check if directory and phonon pickle is here', ...
    'target = this.UserData.dir;', ...
    'if ~isdir(target), target = tempname; mkdir(target); this.UserData.dir=target; end', ...
    'if isempty(dir(fullfile(target, ''phonon.pkl'')))', ...
    '  fid=fopen(fullfile(target, ''phonon.pkl''), ''w'');', ...
    '  if fid==-1, error([ ''model '' this.Name '' '' this.Tag '' could not write phonon.pkl into '' target ]); end', ...
    '  fprintf(fid, ''%s\n'', this.UserData.phonons);', ...
    '  fclose(fid);', ...
    'end', ...
    '  fid=fopen(fullfile(target,''sqw_phonons_eval.py''),''w'');', ...
  [ '  fprintf(fid, ''# Script file for Python/ASE to compute the modes from ' configuration ' in %s\n'', target);' ], ...
    '  fprintf(fid, ''#   ASE: S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002\n'');', ...
    '  fprintf(fid, ''#   <https://wiki.fysik.dtu.dk/ase>. Built by ifit.mccode.org/Models.html sqw_phonons'');', ...
    '  fprintf(fid, ''from ase.phonons import Phonons\n'');', ...
    '  fprintf(fid, ''from os import chdir\n'');', ...
    '  fprintf(fid, ''import numpy\n'');', ...
    '  fprintf(fid, ''import pickle\n'');', ...
    '  fprintf(fid, ''import scipy.io as sio\n'');', ..., ...
    '  fprintf(fid, [ ''chdir("'' target ''")\n'' ]);', ...
    '  fprintf(fid, ''# restore Phonon model\n'');', ...
    '  fprintf(fid, ''fid = open("phonon.pkl", "rb")\n'');', ...
    '  fprintf(fid, ''ph = pickle.load(fid)\n'');', ...
    '  fprintf(fid, ''fid.close()\n'');', ...
    '  fprintf(fid, ''# read HKL locations\n'');', ...
    '  fprintf(fid, ''HKL = numpy.loadtxt("HKL.txt")\n'');', ...
    '  fprintf(fid, ''# compute the spectrum and eigenvectors (polarisation)\n'');', ...
    '  fprintf(fid, ''omega_kn = ph.band_structure(HKL)\n'');', ...
    '  fprintf(fid, ''omega_kn *= 1000\n'');', ...
    '  fprintf(fid, ''# save the result in FREQ\n'');', ...
    '  fprintf(fid, ''sio.savemat("FREQ.mat", { "FREQ": omega_kn })\n'');', ...
    '  fprintf(fid, ''# numpy.savetxt("POLAR", polar_kn)\n'');', ...
    'if isfield(this.UserData, ''DOS'') && isempty(this.UserData.DOS)', ...
    '  fprintf(fid, ''# Calculate phonon DOS\n'');', ...
    '  fprintf(fid, ''print "Computing phonon DOS\\n"\n'');', ...
    '  fprintf(fid, ''omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=5e-4)\n'');', ...
    '  fprintf(fid, ''omega_e *= 1000\n'');', ...
    '  fprintf(fid, ''# save the result in DOS\n'');', ...
    '  fprintf(fid, ''sio.savemat("DOS.mat", { "energy":omega_e, "DOS":dos_e })\n'');', ...
    'end', ...
    '  fprintf(fid, ''exit()\n'');', ...
    '  fclose(fid);', ...
    '  sz0 = size(t);', ...
    '  if ndims(x) == 4, x=squeeze(x(:,:,:,1)); y=squeeze(y(:,:,:,1)); z=squeeze(z(:,:,:,1)); t=squeeze(t(1,1,1,:)); end',...
    'try', ...
    '  if all(cellfun(@isvector, {x y z t})) && numel(unique(cellfun(@(x)length(x), {x y z t}))) > 1',...
    '    [x,y,z] = ndgrid(x,y,z); sz0=[ size(x) sz0(4) ]; end', ...
    '  if all(cellfun(@isscalar,{x y z t})), HKL = [ x y z ; x y z ];', ...
    '  else HKL = [ x(:) y(:) z(:) ]; end', ...
    '  save(''-ascii'', fullfile(target,''HKL.txt''), ''HKL'');', ...
  [ '  [status,result] = system([ ''' precmd 'python '' fullfile(target, ''sqw_phonons_eval.py'') ]);' ], ...
    '  if numel(result) < 1e3, disp(result); end', ...
    '  neg_w = strfind(result,''imaginary frequencies''));', ...
    '  if ~isempty(neg_w), disp([ ''WARNING: found '' num2str(neg_w) '' negative/imaginary phonon frequencies ('' num2str(neg_w/size(HKL,1)) ''%% of total)'' ]); end', ...
    '  clear HKL', ...
    '  delete(fullfile(target,''HKL.txt''))', ...
    '  % import FREQ', ...
    '  FREQ=load(fullfile(target,''FREQ.mat'')); % in meV', ...
    '  FREQ=FREQ.FREQ;', ...
    'catch ME; disp([ ''Model '' this.Name '' '' this.Tag '' could not run Python/ASE from '' target ]); disp(getReport(ME))', ...
    'end', ...
    'try', ...
    'if isfield(this.UserData, ''DOS'') && isempty(this.UserData.DOS)', ...
    '  DOS   = load(fullfile(target,''DOS''),''-ascii'');', ...
    '  DOS_w = load(fullfile(target,''DOS_w''),''-ascii''); DOS=iData(DOS_w,DOS./sum(DOS));', ...
    '  DOS.Title = [ ''DOS '' strtok(this.Name) ]; xlabel(DOS,''DOS Energy [meV]'');', ...
    '  DOS.Error=0; this.UserData.DOS=DOS;', ...
    'end', ...
    'catch ME; disp([ ''Model '' this.Name '' '' this.Tag '' could not get vDOS Python/ASE from '' target ]); disp(getReport(ME));', ...
    'end', ...
    '  % multiply all frequencies(columns, meV) by a DHO/meV', ...
    '  Amplitude = p(1); Gamma=p(2); Bkg = p(3); T=p(4);', ...
    '% apply DHO on each',...
    '  if T<=0, T=300; end', ...
    '  if all(cellfun(@isscalar,{x y z t})), FREQ=FREQ(1,:); end', ...
    '  if all(cellfun(@isvector,{x y z t}))', ...
    '  w=t(:); else w=t(:) * ones(1,size(FREQ,1)); end', ...
    '  signal=zeros(size(w));', ...
    '  this.UserData.maxFreq=max(FREQ(:));', ...
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

% when model is successfully built, display citations
disp(' ');
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '<a href="matlab:doc(''' mfilename ''')">' mfilename '</a>: Model ' configuration ' built using ' options.calculator ])
else
  disp([ mfilename ': Model ' configuration ' built using ' options.calculator ' [' datestr(now) ']' ])
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
case {'QUANTUMESPRESSO','QE'}
cite{end+1} = ' * PHON:   D. Alfe, Computer Physics Communications 180,2622-2633 (2009).';
cite{end+1} = ' * Quantum Espresso: P. Giannozzi, et al J.Phys.:Condens.Matter, 21, 395502 (2009).';
case {'QUANTUMESPRESSO_ASE','QE_ASE'}
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

% sqw_phonons_htmlreport(fullfile(options.target, 'sqw_phonons.html'), 'done', options, cite, signal);

% handle autoplot option
if options.autoplot
  [f, signal] = sqw_phonons_plot(signal);
  if ishandle(options.gui), delete(options.gui); end
else
  f = [];
end

%sqw_phonons_htmlreport(fullfile(options.target, 'sqw_phonons.html'), 'plot', options, f, signal);
%sqw_phonons_htmlreport(fullfile(options.target, 'sqw_phonons.html'), 'end', options, f, signal);

% ------------------------------------------------------------------------------
function [f, signal] = sqw_phonons_plot(signal)
  if isempty(signal) || ~isfield(signal.UserData, 'options') , return; end
  
  options = signal.UserData.options;
  disp([ mfilename ': Model ' options.configuration ' plotting phonons.' ])
  qh=linspace(0.01,1.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);
  f=iData(signal,[],qh,qk,ql,w);
  g=log(f(1,:, :,:)); 
  
  fig1=figure; 
  plot3(g); axis tight;
  view([38 26]);
  slice(g);
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


function signal=sqw_phonons(configuration, varargin)
% model=sqw_phonons(configuration, calculator, ..., options)
%
%   iFunc/sqw_phonons: computes phonon dispersions using the ASE.
%   A model which computes phonon dispersions from the Hellmann-Feynman forces acting 
%     between atoms. The input argument is any configuration file describing the
%     material, e.g. CIF, PDB, POSCAR, ... supported by ASE.
%   This models can compute the coherent inelastic phonons dispersions for
%     any crystalline (powder or single crystal) material, in the harmonic and
%     adiabatic approxiimation.
%   The phonon spectra is computed using one of the calculator supported by the
%   Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase>.
%
%   Supported calculators are:
%     ABINIT    Plane-wave pseudopotential code
%     ELK       Full Potential LAPW code
%     EMT       Effective Medium Theory calculator (Al,Cu,Ag,Au,Ni,Pd,Pt,H,C,N,O)
%     GPAW      Real-space/plane-wave/LCAO PAW code
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
% The model must first be created (which triggers e.g. a DFT computation), and 
% once created, it can be evaluated in the whole HKLE reciprocal space.
%
% MODEL CREATION ===============================================================
%
% sqw_phonons(configuration, calculator, options, ...) 
%   is the usual syntax, such as sqw_phonons('POSCAR','qe','metal')
% sqw_phonons('gui')
%   shows a dialogue to select the phonon calculator configuration.
% sqw_phonons('defaults') or sqw_phonons
%   uses a simple fcc Aluminium as example.
% sqw_phonons(dir, calculator, options, ...)
%   re-use an existing directory containing any of POSCAR, ASE atom object, PhonoPy
%   files, ...
%
% The arguments for the model creation should be:
%
% input (model creation):
% configuration: file name or directory to an existing material configuration
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
% 'html' or 'report':   generate a full report of the simulation and export data.
% 'plot' or 'autoplot': plot the results after building the model.
%
% options: an optional structure with optional settings, as follows:
%
% General options
%   options.target =path                   Where to store all files and FORCES
%     a temporary directory is created when not specified.
%   options.supercell=scalar or [nx ny nz] supercell size. Default is 0 (auto mode).
%   options.calculator=string              EMT,GPAW,Elk,ABINIT,Quantum
%     We recommend ABINIT,QE and Elk. Default set from installed software.
%   options.potentials=string              Basis datasets or pseudopotentials.
%     Elk:    path to potentials, e.g. /usr/share/elk-lapw/species (ELK_SPECIES_PATH)
%     ABINIT: path to potentials, e.g. /usr/share/abinit/psp/      (ABINIT_PP_PATH)
%       or any of 'NC' 'fhi' 'hgh' 'hgh.sc' 'hgh.k' 'tm' 'paw'
%     QuantumEspresso: path to potentials, e.g. /usr/share/espresso/pseudo
%   options.command='exe'                  Path to calculator executable.
%   options.mpi=scalar                     use multi-cpu, requires e.g. OpenMPI to be installed.
%     when available and not specified, MPI is used with all CPU's.
%   options.machinefile=filename           file containing the list of MPI machines to use
%   options.autoplot=0|1                   when set, automatically evaluates the Model
%     and plot results.
%   options.htmlreport=0|1                 when set, automatically generates a full
%     report on the computation results.
%   options.optimizer                      when set, performs a geometry optimization
%     before computing the forces (move atoms in the cell). Can be set to 
%       'BFGS' or 'LBFGS' (low memory BFGS)
%       'MDMin' or 'FIRE'
%   options.accuracy                       'fast', 'very fast' (default) or 'accurate'.
%     The 'fast' choice uses the symmetry operators to lower the number of atom 
%     displacements. The force gradient uses central difference.
%     The 'very fast' option halves the number of displacements. The equilibrium
%     forces are used to improve the accuracy of the force gradient.
%     The 'accurate' choice is longer to execute (e.g. 3-6 times slower), but computes all forces.
%   options.use_phonopy=0|1                requests to use PhonoPy when installed.
%     This choice also sets accuracy='very fast'
%   options.disp=value                     the atom displacement in Angs. Default is 0.01.
%
% DFT specific options
%   options.kpoints=scalar or [nx ny nz]   Monkhorst-Pack grid. Default is 0 (auto mode).
%   options.xc=string                      Type of Exchange-Correlation functional to use
%     'LDA','PBE','revPBE','RPBE','PBE0','B3LYP'   for GPAW
%     'LDA','PBE','revPBE','RPBE','GGA'            for ABINIT
%     'LDA','PBE','REVPBE','PBESOL','WC06','AM05'  for ELK
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
%     selfconsistency, e.g. 1e-8 [eV].
%
% Options specific per calculator
%   options.mode='pw','fd', or 'lcao'      GPAW computation mode as Plane-Wave,
%     Finite Difference, or LCAO (linear combination of atomic orbitals). Default is 'pw'.
%   options.iscf='NC','PAW'                Type of SCF cycles (ABINIT) 
%   options.pps = 'fhi' 'hgh' 'hgh.sc' 'hgh.k' 'tm' 'paw' 'pawxml' Type of database (ABINIT)
%                 a structure such as struct('Ca','_pv') ... (VASP)
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
% The syntax sqw_phonons(model,'html') allows to re-create the HTML report about
% the 4D phonon model.
%
% You may look at the following resources to get material structure files:
%   <http://phonondb.mtl.kyoto-u.ac.jp/>
%   <https://www.materialsproject.org/>
%   <http://nomad-repository.eu/cms/>
%
% WARNING: 
%   This model is suitable to compute phonon dispersions for e.g solid-
%   state materials.
%   The Atomic Simulation Environment must be installed.
%   The temporary directories (UserData.dir) are not removed.
%   Momentum is expressed in reduced lattice unit, and energy is in meV.
%        1 meV = 241.8 GHz = 11.604 K = 0.0965 kJ/mol = 8.0657 cm-1 
%
% References: https://en.wikipedia.org/wiki/Phonon
%
% Atomic Simulation Environment
%   S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002
%   <https://wiki.fysik.dtu.dk/ase>. Exists as Debian package 'python-ase'.
%   Repository: http://download.opensuse.org/repositories/home:/dtufys/
% GPAW J. J. Mortensen et al, Phys Rev B, Vol. 71, 035109 (2005). 
%   <http://wiki.fysik.dtu.dk/gpaw>. Exists as Debian package 'gpaw' and 'gpaw-setups' (gpaw-data).
% Elk <http://elk.sourceforge.net>. Exists as Debian package 'elk-lapw'. 
%   The Elk executable should be compiled with e.g. 
%     elk/src/modmain.f90:289 maxsymcrys=1024 or larger
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
% PhonoPy:   A. Togo and I. Tanaka, Scr. Mater., 108, 1-5 (2015)
%   <https://atztogo.github.io/phonopy/>
%
% MODEL EVALUATION (once created) ==============================================
%
% Once the model has been created, its use requires that axes are given on
%   regular qx,qy,qz grids (in rlu along reciprocal axes). The model evaluations
%   does not require to recompute the forces, and is very fast. 
% When all axes are vectors of same orientation, the HKL locations is assumed to be a q-path.
% When axes are not all vectors, not same length, nor orientation, a 3D HKL cube is used.
% To generate a powder 2D S(q,w) you may use: sqw_powder(model)
% To generate the dispersion curves along principal crystallographic directions,
%   use: sqw_kpath(model)
%
% The phonon density of states (DOS) is automatically computed when evaluating the 
% model onto a grid HKL=[-.5 : 0.5]. It is stored in model.UserData.DOS as an 
% iData object. Partial density of states (per mode) are stored in UserData.DOS_partials
%
% The neutron scattering cross sections should be automatically determined from 
% the chemical formula in order to properly evaluate the neutron scattering 
% intensity. 
% In case this would fail, or other cross sections are needed, it is required
% to set the coherent neutron scattering length (b_coh in [fm]) as a vector in the 
%   model.UserData.properties.b_coh = [vector / number of atoms in the cell]
% or alternatively the scattering cross section in [barn]
%   model.UserData.properties.sigma_coh = [vector / number of atoms in the cell]
% where negative values set b_coh < 0 (e.g. Hydrogen).
% Setting b_coh=0 unactivates the intensity estimate.
% This assignement should be done after creating the model, before performing
% further HKL evaluations. Default is to use guessed b_coh from the formula.
%
% Once the model has been created:
% input:  p: sqw_phonons model parameters (double)
%             p(1)=Amplitude
%             p(2)=Gamma   dispersion DHO half-width in energy [meV]
%             p(3)=Background (constant)
%             p(4)=Temperature of the material [K]. When 0, the intensity is not computed.
%             p(5)=Energy_scaling. All frequencies are multiplied by this value.
%          or p='guess'
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
%
% Example (model creation and evaluation):
%   s=sqw_phonons('bulk("Cu", "fcc", a=3.6, cubic=True)','EMT','metal');
%   qh=linspace(-.5,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);
%   f=iData(s,[],qh,qk,ql',w); scatter3(log(f(1,:, :,:)),'filled');
%   figure; plot(s.UserData.DOS); % plot the DOS
%
%   sqw_phonons(s, 'report') generates a full report for the given model.
%   sqw_kpath(s)             plots the dispersion curves and density of states.
%
%   s=sqw_phonons('bulk("Si", "diamond", a=5.4, cubic=True)','semiconductor');
%
%   s=sqw_phonons([ ifitpath 'Data/POSCAR_Al'],'metal','EMT');
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_cubic_monoatomic, sqw_sine3d, sqw_vaks
%   sqw_powder, <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

% Units: 1 Ry        = 13.6 eV
%        1 Ha = 2 Ry = 27.2 eV

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
options.available = status;
% check if we re-use an existing iFunc Model
if isa(configuration, 'iFunc') && configuration.Dimension == 4
  signal = configuration;
  if ~options.htmlreport, return; end

  if isfield(configuration.UserData,'options')
    options = configuration.UserData.options;
  end
  options.htmlreport = 1;
  % cope with other phonon models
  if ~isfield(options, 'configuration')
    options.configuration = signal.Description;
  end
  if ~isfield(options, 'calculator')
    options.calculator = signal.Name;
  end
  if ~isfield(options, 'duration')
    options.duration = 0;
  end
  Phonon_Model = signal;
  if ~isdir(options.target), mkdir(options.target); end
  builtin('save', fullfile(options.target, 'Phonon_Model.mat'), 'Phonon_Model');

  sqw_phonons_htmlreport('', 'create_atoms', options);
  sqw_phonons_htmlreport('', 'results', options);
  filename = sqw_phonons_htmlreport('', 'download', options);

  signal.UserData.options.report = filename;
  signal.UserData.options.target = options.target;
  if ~isempty(inputname(1))
    assignin('caller',inputname(1),signal);
  end
  return
end

% check if we point to an already existing directory with previous Model/FORCES
if ~isempty(dir(configuration)) % a file/directory
  % try to load the configuration as a Mat file
  try
    signal = load(configuration);
  end
  if isempty(signal) && isdir(configuration)
    % look for Phonon_Model.mat
    if ~isempty(dir(fullfile(configuration,'Phonon_Model.mat')))
      signal = load(fullfile(configuration,'Phonon_Model.mat'));
    end
  end
  
  if isa(signal, 'iFunc') && ndims(signal) == 4
    return
  end
  if isstruct(signal)
    % search for a 4D model, and return if found
    for f=fieldnames(signal)
      this = signal.(f{1});
      if isa(this, 'iFunc') && ndims(this) == 4
        signal = this;
        return
      end
    end
  end
end

if isempty(status.mpirun) && isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
  options.mpi=1;
  disp([ mfilename ': MPI parallelization is not available. Install e.g. OpenMPI first. Using mpi=1 (serial).' ]);
end

if isempty(options.calculator)
  % select default calculator according to those installed
  for index={'quantumespresso','quantumespresso_ase',  ...
      'abinit','vasp','gpaw','elk'};
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
  options = sqw_phonons_gui(configuration, options, status);
  if isempty(options), return; end
end % GUI

% make sure we find the 'configuration' also from ifitpath
if isempty(dir(configuration)) && ~isempty(dir(fullfile(ifitpath,'Data',configuration)))
  configuration = fullfile(ifitpath,'Data',configuration);
end

% ==============================================================================
%                             BUILD MODEL (read initial structure)
%                                    sqw_phonons_check (.py)
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

% BUILD stage: we call ASE to build the model
target = options.target;
[options, result, read] = sqw_phonons_check(configuration, options, status);
configuration = options.configuration; % in case this has been updated in 'check'

% we test if the pickle file could be written. This way even if the export/save 
% properties fail, we can proceed.
if isempty(dir(fullfile(target, 'atoms.pkl')))  % FATAL
  sqw_phonons_error([ mfilename ': failed read input ' configuration ], options);
  return
end

options.report = sqw_phonons_htmlreport('', 'create_atoms', options);

% ==============================================================================
%                               BUILD MODEL (get calculator)
%                                     sqw_phonons_calc
% ==============================================================================

% get the calculator
% for QE, this triggers computation with sqw_phon()
[decl, calc, signal,options] = sqw_phonons_calc(options, status, options.calculator, read);
% return directly is the signal=iFunc has been created (QE case)
if isempty(decl) && isempty(signal), return; end

% ==============================================================================
%                               BUILD MODEL (optimize or pass)
%                                   sqw_phonons_optimize (.py)
% ==============================================================================

% handle the optimizer
% requires 'atoms' and 'calc'. provides new 'atoms'.

[options, sav] = sqw_phonons_optimize(decl, calc, options);

% ==============================================================================
%                               COMPUTE MODEL (force computation)
%                                   sqw_phonons_get_forces
%                               sqw_phonons_forces_iterate.py
%                               sqw_phonons_forces_finalize.py
%                                    not for QE/PHON case
% ==============================================================================

% the QE case has been done prior to optimize when calling sqw_phonons_calc()

if ~strcmpi(options.calculator, 'QUANTUMESPRESSO') || strcmpi(options.calculator, 'QUANTUMESPRESSO_ASE')

  [options, sav] = sqw_phonons_get_forces(options, decl, calc);
  if isempty(options), return; end

  % create the iFunc object
  % then read the pickle file to store it into the model
  if ~isempty(dir(fullfile(target, 'phonon.pkl')))
    signal.UserData.phonons = binary_fileread(fullfile(target, 'phonon.pkl')); % binary
  else
    if ~isdeployed && usejava('jvm') && usejava('desktop')
      target = [ '<a href="' target '">' target '</a>' ];
    end
    sqw_phonons_error(sprintf([ mfilename ': ' options.calculator ...
      ' failed. Temporary files and Log are in ' target ...
      '\n  Try increasing the number of iterations nsteps=' num2str(max(200, options.nsteps*2)) ...
      '\n  or a relaxed tolerance toldfe=' num2str(max(1e-5, options.toldfe*10)) ...
      '\n  or there is a missing pseudopotential, or other stuff...' ]), options)
    sqw_phonons_htmlreport('', 'download', options);
    return
  end
  if ~isempty(dir(fullfile(target, 'phonopy.pkl')))
    signal.UserData.phonopy = binary_fileread(fullfile(target, 'phonopy.pkl')); % binary
  end
  
  if ~isempty(dir(configuration))
    [dummy, signal.UserData.input]= fileparts(configuration);
  else
    signal.UserData.input = configuration;
  end

  signal.Parameters     = {  ...
    'Amplitude' ...
    'Gamma Damped Harmonic Oscillator width in energy [meV]' ...
    'Background' ...
    'Temperature used to compute the Bose factor n(w) [K]' ...
    'Energy_scaling' ...
     };
    
  signal.Dimension      = 4;         % dimensionality of input space (axes) and result

  signal.Guess = [ 1 .1 0 10 1 ];
  
  if ~isempty(dir(fullfile(target, 'FORCE_SETS')))  % always prefer to get PhonoPy when exists
    signal.UserData.FORCES = fileread(fullfile(target, 'FORCE_SETS'));
  elseif ~isempty(dir(fullfile(target, 'FORCES.mat')))
    signal.UserData.FORCES = load(fullfile(target, 'FORCES.mat'));
  end
  
  % get 'atoms' back from python and compute neutron cross sections from formulae
  if ~isempty(fullfile(target, 'properties.mat'))
    properties = load(fullfile(target, 'properties.mat'));
    if isfield(properties, 'chemical_symbols')
      properties.b_coh = sqw_phonons_b_coh(properties.chemical_symbols);
      save(fullfile(target, 'properties.mat'), '-v6', 'properties')
    end
    if isfield(properties, 'masses')
      properties.molar_mass = sum(properties.masses);
    end
    
    % compute the density
    if isfield(properties, 'masses') && isfield(properties, 'volume')
      mole = 6.022140857e23;
      properties.density = sum(properties.masses)/properties.volume/mole*1e24;
    end

  else
    properties = [];
  end
  
  if ~isempty(dir(configuration)) && ~isdir(configuration)
    signal.UserData.configuration = fileread(configuration);
  else
    signal.UserData.configuration = configuration;
  end

  signal.UserData.dir           = target;
  signal.UserData.options       = options;
  signal.UserData.calc          = calc;
  signal.UserData.properties    = orderfields(properties);
  
  if options.use_phonopy
    c = [ 'PhonoPy/' options.calculator ];
  else
    c = [ 'ASE/' options.calculator ];
  end
  if isfield(properties, 'chemical_formula'), m=properties.chemical_formula;
  else                                        m=signal.UserData.input; end
  
  
  signal.Description    = [ 'S(q,w) ' m ' dispersion(HKL) Phonon/' c ' with DHO line shape. ' configuration ];
  signal.Name           = [ 'Sqw_Phonon_' m ' ' c ' DHO [' mfilename ']' ];

  % EVAL stage: we call ASE to build the model. ASE does not support single HKL location. 
  % For this case we duplicate xyz and then get only the 1st FREQ line
  
  script_restore_pickles = { ...
      'if isempty(dir(fullfile(target, ''phonon.pkl'')))', ...
    '  fid=fopen(fullfile(target, ''phonon.pkl''), ''w'');', ...
    '  if fid==-1, error([ ''model '' this.Name '' '' this.Tag '' could not write phonon.pkl into '' target ]); end', ...
    '  fwrite(fid,this.UserData.phonons);', ...
    '  fclose(fid);', ...
    'end', ...
    'if isempty(dir(fullfile(target, ''phonopy.pkl''))) && isfield(this.UserData,''phonopy'')', ...
    '  fid=fopen(fullfile(target, ''phonopy.pkl''), ''w'');', ...
    '  if fid==-1, error([ ''model '' this.Name '' '' this.Tag '' could not write phonopy.pkl into '' target ]); end', ...
    '  fwrite(fid, this.UserData.phonopy);', ...
    '  fclose(fid);', ...
    'end', ...
    'if isempty(dir(fullfile(target, ''ifit.py'')))', ...
    '  fid=fopen(fullfile(target, ''ifit.py''), ''w'');', ...
    '  if fid==-1, error([ ''model '' this.Name '' '' this.Tag '' could not write ifit.py into '' target ]); end', ...
    '  fprintf(fid, ''%s\n'', this.UserData.options.script_ifitpy);', ...
    '  fclose(fid);', ...
    'end', ...
  };
  
  % we use the template 'snippets' which are in the private directory.
  
  % get code to read xyzt and build HKL list and convolve DHO line shapes
  [script_hkl, script_dho] = sqw_phonons_templates;
  
  script_compute_HKL = { ...
      '  if ~this.UserData.options.use_phonopy || isempty(this.UserData.options.available.phonopy)', ...
      '    fprintf(fid, ''fid = open("phonon.pkl", "rb")\n'');', ...
      '    fprintf(fid, ''ph = pickle.load(fid)\n'');', ...
      '    fprintf(fid, ''fid.close()\n'');', ...
      '    fprintf(fid, ''omega_kn, polar_kn = ph.band_structure(HKL, modes=True,verbose=False)\n'');', ...
      '    fprintf(fid, ''omega_kn *= 1000\n'');', ...
      '  else', ...
      '    fprintf(fid, ''fid = open("phonopy.pkl", "rb")\n'');', ...
      '    fprintf(fid, ''ph = pickle.load(fid)\n'');', ...
      '    fprintf(fid, ''fid.close()\n'');', ...
      '    fprintf(fid, ''omega_kn = []\n'');', ...
      '    fprintf(fid, ''omega_kn, polar_kn = ifit.phonopy_band_structure(ph, HKL, modes=True)\n'');', ...
      '    fprintf(fid, ''omega_kn *= 4.13 # from THz to meV\n'');', ...
      '  end' ...
    };
  
  % python scripts for band structure evaluation
  % initially we saved using:
  %   sio.savemat("FREQ.mat", { "FREQ": omega_kn, "POLAR":polar_kn, "HKL":HKL })
  % but these MAT files are limited to 2Gb.
  % We now use hdf5storage.savemat
  signal.Expression = { ...
    '% check if directory and phonon pickle is here', ...
    'target = this.UserData.dir;', ...
    'if ~isdir(target), target = tempname; mkdir(target); this.UserData.dir=target; end', ...
    script_restore_pickles{:}, ...
    '  fid=fopen(fullfile(target,''sqw_phonons_eval.py''),''w'');', ...
  [ '  fprintf(fid, ''# Script file for Python/ASE to compute the modes from ' configuration ' in %s\n'', target);' ], ...
    '  fprintf(fid, ''#   ASE: S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002\n'');', ...
    '  fprintf(fid, ''#   <https://wiki.fysik.dtu.dk/ase>. Built by ifit.mccode.org/Models.html sqw_phonons'');', ...
    '  fprintf(fid, ''from ase.phonons import Phonons\n'');', ...
    '  fprintf(fid, ''from os import chdir\n'');', ...
    '  fprintf(fid, ''import numpy\n'');', ...
    '  fprintf(fid, ''import pickle\n'');', ...
    '  fprintf(fid, ''import ifit\n'');', ...
    '  fprintf(fid, ''import hdf5storage\n'');', ...
    '  fprintf(fid, [ ''chdir("'' target ''")\n'' ]);', ...
    '  fprintf(fid, ''# restore Phonon model and read HKL locations\n'');', ...
    '  fprintf(fid, ''HKL = numpy.loadtxt("HKL.txt")\n'');', ...
    '  fprintf(fid, ''# compute the spectrum and eigenvectors (polarisation)\n'');', ...
    script_compute_HKL{:}, ...
    '  fprintf(fid, ''# save the result in FREQ into MAT/HDF5 file\n'');', ...
    '  fprintf(fid, ''hdf5storage.savemat("FREQ.mat", { "FREQ":omega_kn, "POLAR":polar_kn, "HKL":HKL })\n'');', ...
    '  fprintf(fid, ''exit()\n'');', ...
    '  fclose(fid);', ...
    script_hkl{:}, ...
    '% ASE does not support single HKL location. We duplicate it when found.', ...
    'if size(HKL,1) == 1, flag_single_point=1; HKL = [ HKL ; HKL ]; else flag_single_point=0; end', ...
    'try', ...
    '  save(''-ascii'', fullfile(target,''HKL.txt''), ''HKL'');', ...
  [ '  [status,result] = system([ ''' precmd 'python '' fullfile(target, ''sqw_phonons_eval.py'') ]);' ], ...
    '  if numel(result) > 1e3, disp(result(1:500)); disp(''...''); disp(result((end-480):end)); ', ...
    '  else disp(result); end', ...
    '  clear result', ...
    '  % import FREQ', ...
    '  FREQ=load(fullfile(target,''FREQ.mat''));', ...
    '  delete(fullfile(target,''FREQ.mat''));', ...
    '  POLAR=FREQ.POLAR; FREQ=FREQ.FREQ;', ...
    'catch ME; disp([ ''Model '' this.Name '' '' this.Tag '' could not run Python/ASE from '' target ]); disp(getReport(ME));', ...
    'end', ...
    'if flag_single_point, FREQ=FREQ(1,:); POLAR=POLAR(1,:); HKL=HKL(1,:); end', ...
    '% multiply all frequencies(columns, meV) by a DHO/meV', ...
    'Amplitude = p(1); Gamma=p(2); Bkg = p(3); T=p(4); ', ...
    'if p(5) > 0, FREQ=FREQ*p(5); end', ...
    script_dho{:} };

  signal = iFunc(signal);
end % other calculators than QE

cite = { ' * Atomic Simulation Environment', ...
  '           S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002.', ...
  ' * iFit:   E. Farhi et al, J. Neut. Res., 17 (2013) 5.' };
switch upper(options.calculator)
case 'GPAW'
cite{end+1} = ' * GPAW:   J. J. Mortensen et al, Phys Rev B, Vol. 71, 035109 (2005).';
case 'ELK'
cite{end+1} = ' * ELK:    http://elk.sourceforge.net';
case 'ABINIT'
cite{end+1} = ' * ABINIT: X. Gonze et al, Computer Physics Communications 180, 2582-2615 (2009).';
case 'EMT'
cite{end+1} = ' * EMT:    K.W. Jacobsen et al, Surf. Sci. 366, 394-402 (1996).';
case {'QUANTUMESPRESSO','QE'}
cite{end+1} = ' * Quantum Espresso: P. Giannozzi, et al J.Phys.:Condens.Matter, 21, 395502 (2009).';
cite{end+1} = ' * PHON:   D. Alfe, Computer Physics Communications 180,2622-2633 (2009).';
case {'QUANTUMESPRESSO_ASE','QE_ASE'}
cite{end+1} = ' * Quantum Espresso: P. Giannozzi, et al J.Phys.:Condens.Matter, 21, 395502 (2009).';
cite{end+1} = ' * QE-util:          Pawel T. Jochym, https://jochym.github.io/qe-doc/ (2015).';
case 'VASP'
cite{end+1} = ' * VASP:   G. Kresse and J. Hafner. Phys. Rev. B, 47:558, 1993.';
case 'OCTOPUS'
cite{end+1} = ' * OCTOPUS:X. Andrade et al, Phys Chem Chem Phys 17 (2015) 31371.';
end
if options.use_phonopy
  cite{end+1} = ' * PhonoPy:   A. Togo and I. Tanaka, Scr. Mater., 108, 1-5 (2015)';
end

signal.UserData.duration = etime(clock, t);
options.duration = signal.UserData.duration;
options.cite     = cite;
signal.UserData.options   = orderfields(options);
signal.UserData  = orderfields(signal.UserData);

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

disp([ 'Time elapsed=' num2str(signal.UserData.duration) ' [s]. Please cite:' ])

fprintf(1, '%s\n', cite{:});
disp('You can now evaluate the model using e.g.:')
disp('    qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);');
disp('    f=iData(s,[],qh,qk,ql'',w); plot3(log(f(1,:, :,:)));');
disp(' ');

% save the Model as a Matlab object
Phonon_Model = signal;
builtin('save', fullfile(options.target, 'Phonon_Model.mat'), 'Phonon_Model');

options.status = 'DONE';
sqw_phonons_htmlreport('', 'status', options);

sqw_phonons_htmlreport('', 'results', options);

% handle autoplot option
if options.autoplot
  [S, k, fig] = sqw_kpath(Phonon_Model,'','','plot THz');
  if ishandle(options.gui), delete(options.gui); end
else
  f = [];
end

sqw_phonons_htmlreport('', 'download', options);

% =============================================================================
function content = binary_fileread(filename)
  content = [];
  fid = fopen(filename, 'r');
  if fid == -1, return; end
  content = uint8(fread(fid));
  fclose(fid);

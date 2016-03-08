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
%     NWChem    Gaussian based electronic structure code (with MPI)
%     Dacapo    Plane-wave ultra-soft pseudopotential code
%     ELK       Full Potential LAPW code (with MPI)
%     ABINIT    Plane-wave pseudopotential code (with MPI)
%     QuantumEspresso Plane-wave pseudopotential code (with MPI)
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
% The arguments for the model creation should be:
%
% configuration: file name to an existing material configuration
%   Any A.S.E supported format can be used (POSCAR, CIF, SHELX, PDB, ...). 
%     See <https://wiki.fysik.dtu.dk/ase/ase/io.html#module-ase.io>
%   Alternatively, the 'bulk','molecule', and 'nanotube' ASE constructors can be
%   used, using the Python syntax, e.g. 
%       'bulk("Si", "diamond", a=5.4)'
%       'bulk("Cu", "fcc", a=3.6, cubic=True)'
%       'molecule("H2O")'
%       'nanotube(6, 0, length=4)'
%     See <https://wiki.fysik.dtu.dk/ase/ase/structure.html>
%
% 'metal','insulator','semiconductor': indicates the type of occupation for 
%    electronic states, which sets smearing.
%
% options: an optional structure with optional settings
%
% General options
%   options.target =path                   Where to store all files and FORCES
%     a temporary directory is created when not specified.
%   options.supercell=scalar or [nx ny nz] supercell size. Default is 2.
%   options.calculator=string              EMT,GPAW,Elk,NWChem,Dacapo,ABINIT,Quantum
%     We recommend ABINIT,QE and Elk. Default set from installed software.
%   options.dos=1                          Option to compute the vibrational
%     density of states (vDOS) in UserData.DOS
%   options.potentials=string              Basis datasets or pseudopotentials.
%     NWChem: see http://www.nwchem-sw.org/index.php/Release64:AvailableBasisSets
%       e.g. 'ANO-RCC' 'ANO-RCC' 
%     Dacapo: path to potentials, e.g. /usr/share/dacapo-psp       (DACAPOPATH)
%     Elk:    path to potentials, e.g. /usr/share/elk-lapw/species (ELK_SPECIES_PATH)
%     ABINIT: path to potentials, e.g. /usr/share/abinit/psp/      (ABINIT_PP_PATH)
%       or any of 'NC' 'fhi' 'hgh' 'hgh.sc' 'hgh.k' 'tm' 'paw'
%     QuantumEspresso: path to potentials, e.g. /usr/share/espresso/pseudo
%   options.command='exe'                  Path to calculator executable, which may
%     include e.g. "mpirun -np N" when applicable.
%   options.mpi=scalar                     use multi-cpu, for NWChem,QuantumEspresso
%     which requires e.g. OpenMPI to be installed.
%
% DFT specific options
%   options.kpoints=scalar or [nx ny nz]   Monkhorst-Pack grid
%   options.xc=string                      Type of Exchange-Correlation functional to use
%     'LDA','PBE','revPBE','RPBE','PBE0','B3LYP'   for GPAW
%     'LDA','B3LYP','PBE','RHF','MP2'              for NWChem
%     'LDA','PBE','REVPBE','PBESOL','WC06','AM05'  for ELK
%     'PZ','VWN','PW91','PBE','RPBE','revPBE'      for Dacapo/Jacapo
%     'LDA', 'PBE', 'revPBE', 'RPBE'               for ABINIT
%     for QuantumEspresso, this choice is made through a dialog.
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
%   options.mixing_beta=scalar             mixing factor for self-consistency
%     default=0.7. use 0.3 to improve convergence (QuantumEspresso)
%   options.mixing_ndim=scalar             number of iterations used in mixing
%     default=8. If you are tight with memory, you may reduce it to 4. A larger
%     value will improve the SCF convergence, and use more memory (QuantumEspresso)
%
% The options can also be entered as strings with 'field=value; ...'.
%
% Once the model has been created, its use requires that axes are given on
% regular qx,qy,qz grids.
%     
% Example:
%   s=sqw_phonons('bulk("Cu", "fcc", a=3.6, cubic=True)','EMT','metal');
%   qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);
%   f=iData(s,[],qh,qk,ql,w); scatter3(log(f(1,:, :,:)),'filled');
%   figure; plot(s.UserData.DOS); % plot the DOS, as indicated during model creation
%
%   s=sqw_phonons('bulk("Si", "diamond", a=5.4, cubic=True)','semiconductor');
%
%   s=sqw_phonons([ ifitpath 'Data/POSCAR_Al'],'dos','metal','EMT');
%
% WARNING: Single intensity and line width parameters are used here.
%   This model is suitable to compute phonon dispersions for e.g solid-
%   state materials.
%   The Atomic Simulation Environment must be installed.
%   The temporary directories (UserData.dir) are not removed.
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
%
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
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>

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

options= sqw_phonons_argin(configuration, varargin{:});

if isempty(options.calculator)
  % select default calculator according to those installed
  for index={'quantumespresso','abinit','gpaw','elk','jacapo','nwchem',};
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

t=clock();

if options.gui
  % pop-up a simple dialog requesting for:
  %  * configuration
  %  * calculator
  %  * metal, insulator, semiconductor
  %  * supercell
  %  * kpoints
  % and will select autoplot, 'dos'
  calcs = 'EMT';
  for index={'gpaw','elk','jacapo','nwchem','abinit','quantumespresso'};
    if ~isempty(status.(index{1})), calcs = [ calcs ', ' upper(index{1}) ]; end
  end
  NL = sprintf('\n');
  prompt = { [ '{\bf Atom/molecule/system configuration}' NL 'a CIF/PDB/POSCAR/... name or e.g. bulk("Cu", "fcc", a=3.6, cubic=True),' NL  'molecule("H2O"), or nanotube(6, 0, length=4). ' NL 'Documentation at {\color{blue}http://ifit.mccode.org/Models.html}' ], ...
  [ '{\bf Calculator}' NL 'One of ' calcs ], ...
  [ '{\bf Smearing}' NL 'metal, semiconductor, insulator or left empty. Use e.g. 0.3 eV for conductors, or a small value such as 0.01 to help SCF convergence. You may use "auto" with Elk. ' ], ...
  [ '{\bf Cut-off energy for wave-functions}' NL 'Leave as 0 for the default, or specify a cut-off in eV, e.g. 500 eV for fast estimate, 1500 or 2000 eV for accurate results.' ], ...
  [ '{\bf K-Points}' NL 'Monkhorst-Pack grid which determines the K sampling (3-vector). 4 is the minimum for accurate computations, 6 is best. Use 1 or 2 for testing only (faster).' ], ...
   [ '{\bf Supercell}' NL 'The size of the repeated model = system*supercell (3-vector). Should be larger than k-points.' ], ...
   [ '{\bf Other options}' NL 'Such as mpi, nbands, nsteps, xc (default PBE), toldfe, raw' NL 'example: "mpi=4; nsteps=50"' NL 'Documentation at {\color{blue}http://ifit.mccode.org/Models.html}' ] };
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
end

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
  end
end

Ha = 27.2; Ry=13.6;





% ------------------------------------------------------------------------------
% handle supported calculators: decl and calc
if options.gui && ishandle(options.gui), waitbar(0.05, options.gui, [ mfilename ': configuring ' options.calculator ]); end
if isunix, precmd = 'LD_LIBRARY_PATH= ; '; else precmd=''; end
switch upper(options.calculator)

% ==============================================================================
case 'ABINIT'
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  if isempty(strfind(status.(lower(options.calculator)),'abinis')) && isempty(options.command)
    options.command = status.(lower(options.calculator));
  end
  if isfield(options,'mpi') && ~isempty(options.mpi)
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    if isscalar(options.mpi) 
      options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
    end
  end
  if ~isempty(options.command)
    cmd = options.command;
    if isempty(strfind(cmd, 'PREFIX.files'))
      cmd = [ cmd ' < PREFIX.files > PREFIX.log' ];
    end
    setenv('ASE_ABINIT_COMMAND', cmd);
  end
  if ~isempty(options.potentials)
    if strcmpi(options.potentials,'NC')
      options.potentials='';
      options.iscf=7;
    elseif any(strcmpi(options.potentials, {'fhi', 'hgh', 'hgh.sc', 'hgh.k', 'tm', 'paw'}))
      if strcmpi(options.potentials, 'paw')
        options.iscf=17; % seems best. see https://www.nsc.liu.se/~pla/blog/2014/02/21/deltacodes/
      end
      options.pps = lower(options.potentials);
      options.potentials = '';
    end
  end
  if isunix
    if isempty(options.potentials), options.potentials='/usr/share/abinit/psp/'; end
  end
  if ~isempty(options.potentials) && isdir(options.potentials)
    setenv('ABINIT_PP_PATH', options.potentials);
    d = [ dir(fullfile(options.potentials,'LDA_*')) ; ...
          dir(fullfile(options.potentials,'GGA_*')) ];
    for index=d'
      if index.isdir
        setenv('ABINIT_PP_PATH', ...
          [ getenv('ABINIT_PP_PATH') pathsep fullfile(options.potentials, index.name) ]); 
      end
    end
  end
  % parallelisation: npbands npftt https://www.nsc.liu.se/~pla/blog/2012/04/18/abinitvasp-part2/
  decl = 'from ase.calculators.abinit import Abinit';
  calc = 'calc = Abinit(chksymbreak=0 ';
  if options.ecut <= 0, options.ecut=340; end % no default in ABINIT (eV)
  if (options.ecut > 0)
    calc = [ calc sprintf(', ecut=%g', options.ecut) ];
  end
  if options.toldfe <= 0, options.toldfe=1e-5; end % in eV, necessary
  if (options.toldfe > 0)
    calc = [ calc sprintf(', toldfe=%g', options.toldfe) ];
  end
  if isfield(options,'iscf')
    % iscf=7 default (NC), 17 (PAW) 
    % <http://www.abinit.org/doc/helpfiles/for-v7.0/input_variables/varbas.html#iscf>
    if options.iscf < 0, options.iscf=7; end 
    if (options.iscf > 0)
      calc = [ calc sprintf(', iscf=%i', options.iscf) ];
    end
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=[%i,%i,%i]', options.kpoints) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if isfield(options,'mpi') && options.mpi > 1
    % nbdblock, npband, AUTOPARAL=1
    % calc = [ calc sprintf(', nbdblock=%i', options.mpi) ];
    calc = [ calc sprintf(', autoparal=1') ];
  end
  if ~isfield(options, 'pps') || isempty(options.pps)
    options.iscf=17;
  end
  if isfield(options, 'pps') && ~isempty(options.pps)
    calc = [ calc sprintf(', pps=''%s''', options.pps) ];
  end
  if options.nbands > 0
    calc = [ calc sprintf(', nband=%i', options.nbands) ];
  end
  if options.nsteps > 0
    calc = [ calc sprintf(', nstep=%i', options.nsteps) ];
  end
  if isscalar(options.occupations) && options.occupations >=0
    calc=[ calc sprintf(', tsmear=%g, occopt=3', options.occupations/Ha) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];

case 'ELK' % ===================================================================
  % requires custom compilation with elk/src/modmain.f90:289 maxsymcrys=1024
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  % location of ELF pseudo-potentials is mandatory
  if isempty(options.potentials) && isempty(getenv('ELK_SPECIES_PATH'))
    if isunix, options.potentials = '/usr/share/elk-lapw/species';
      disp([ mfilename ': ' options.calculator ': assuming atom species are in' ])
      disp([ '  ' options.potentials ])
      disp('  WARNING: if this is not the right location, use options.potentials=<location>');
    else
      sqw_phonons_error([ mfilename ': ' options.calculator ': undefined "species". Use options.potentials=<location of elk/species>.' ], options)
    end
  end
  if ~isempty(options.potentials)
    setenv('ELK_SPECIES_PATH', [ options.potentials, filesep ]);
  end
  if ~strcmp(status.(lower(options.calculator)),'elk') && isempty(options.command)
    options.command = status.(lower(options.calculator));
  end
  if isfield(options,'mpi') && ~isempty(options.mpi)
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    if isscalar(options.mpi) 
      options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
    end
  end
  if ~isempty(options.command)
    cmd = options.command;
    if isempty(strfind(cmd, 'elk.out'))
      cmd = [ cmd ' > elk.out' ];
    end
    setenv('ASE_ELK_COMMAND', cmd);
  end
  
  decl = 'from ase.calculators.elk import ELK';
  calc = 'calc = ELK(tforce=True, tasks=0, rgkmax=4.0, mixtype=3'; % Pulay mixing

  if strcmp(options.occupations, 'auto')
    % other distribution: MethfesselPaxton
    calc=[ calc sprintf(', stype=3, autoswidth=True') ];
  elseif isscalar(options.occupations) && options.occupations >=0
    calc=[ calc sprintf(', stype=3, swidth=%g', options.occupations/Ha) ];
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=(%i,%i,%i)', options.kpoints) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if options.nbands > 0
    calc = [ calc sprintf(', nvbse=%i', options.nbands) ];
  end
  if options.nsteps > 0
    calc = [ calc sprintf(', maxscl=%i', options.nsteps) ];
  end
  if options.toldfe > 0
    calc = [ calc sprintf(', epsengy=%g', options.toldfe) ];
  end
  if options.ecut > 0
    calc = [ calc sprintf(', emaxrf=%g', options.ecut) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  
case 'EMT'
  decl = 'from ase.calculators.emt import EMT';
  calc = 'calc  = EMT()';
case 'GPAW' % ==================================================================
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  
  decl = 'from gpaw import GPAW, PW, FermiDirac';
  calc = 'calc = GPAW(usesymm=False'; % because small displacement breaks symmetry

  if isscalar(options.occupations) && options.occupations>=0 % smearing in eV
    calc=[ calc sprintf(', occupations=FermiDirac(%g)', options.occupations) ];
    % other distribution: MethfesselPaxton
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=(%i,%i,%i)', options.kpoints) ];
  end
  if options.ecut > 0
    calc = [ calc sprintf(', mode=PW(%g)', options.ecut) ];
  elseif ~isempty(options.mode)
    calc = [ calc sprintf(', mode=''%s''', options.mode) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if ~isempty(options.diagonalization)
    if strncmpi(options.diagonalization, 'dav', 3) options.diagonalization='dav'; end
    calc = [ calc sprintf(', eigensolver=''%s''', options.diagonalization) ];
  end
  if ~isempty(options.potentials)
    calc = [ calc sprintf(', setups=''%s''', options.potentials) ];
  end
  if options.nbands > 0
    calc = [ calc sprintf(', nbands=%i', options.nbands) ];
  end
  if options.nsteps > 0
    calc = [ calc sprintf(', maxiter=%i', options.nsteps) ];
  end
  if options.toldfe > 0
    calc = [ calc sprintf(', convergence={''energy'':%g}', options.toldfe) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  
case 'JACAPO' % ================================================================
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  
  % Requires to define variables under Ubuntu
  if isunix
    if isempty(options.potentials), options.potentials='/usr/share/dacapo-psp'; end
    if isempty(options.command),    options.command   =status.(lower(options.calculator)); end
  end
  if ~isempty(options.potentials)
    setenv('DACAPOPATH', options.potentials);
  end
  if isfield(options,'mpi') && ~isempty(options.mpi)
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    if isscalar(options.mpi) 
      options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
    end
  end
  if ~isempty(options.command)
    setenv('DACAPOEXE_SERIAL', options.command); % DACAPOEXE_PARALLEL
  end
  
  decl = 'from ase.calculators.jacapo import Jacapo';
  calc = 'calc = Jacapo(symmetry=False';
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if options.ecut < 0
    % options.ecut = 340;
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=(%i,%i,%i)', options.kpoints) ];
  end
  if options.ecut > 0
    calc = [ calc sprintf(', pw=%g', options.ecut) ];
  end
  if options.nbands > 0
    calc = [ calc sprintf(', nbands=%i', options.nbands) ];
  end
  if options.occupations > 0
    calc = [ calc sprintf(', ft=%g', options.occupations) ];
  end
  if isfield(options,'ecutrho') && optins.ecutrho > 0
    calc = [ calc sprintf(', dw=%g', options.ecutrho) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];  
  % other options
  if optins.toldfe > 0
    calc = [ calc sprintf('; calc.SetConvergenceParameters(%g) ', options.toldfe) ];
  end
  if ~isempty(options.diagonalization)
    if strncmp(options.diagonalization, 'dav',3), options.diagonalization='eigsolve';
    elseif strcmp(options.diagonalization, 'cg'), options.diagonalization='resmin';
    else options.diagonalization='rmm-diis'; end
    calc = [ calc sprintf('; calc.SetEigenvalueSolver(%s)', ...
      options.diagonalization) ];
  end

case 'NWCHEM' % ================================================================
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  if isfield(options,'mpi') && ~isempty(options.mpi)
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    if isscalar(options.mpi) 
      options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
    end
  end
  if ~isempty(options.command)
    cmd = options.command;
    if isempty(strfind(cmd, 'PREFIX.nw'))
      cmd = [ cmd ' PREFIX.nw > PREFIX.out' ];
    end
    setenv('ASE_NWCHEM_COMMAND', cmd);
  end
  
  
  
  decl = 'from ase.calculators.nwchem import NWChem';
  calc = sprintf('calc = NWChem(xc=''%s''', ...
    options.xc);
  % check if we use KPTS
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', raw="nwpw\\n  monkhorst-pack %i %i %i\\nend"', options.kpoints) ];
  end
  if ~isempty(options.potentials)
    calc = [ calc sprintf(', basis=''%s''', options.potentials) ];
  end
  % smearing is in Hartree
  if strcmp(options.occupations, 'smearing') || strcmp(options.occupations, 'metal') % metals
    options.occupations=0.1;
  elseif strcmp(options.occupations, 'semiconductor')
    options.occupations=0;
  elseif strcmp(options.occupations, 'fixed') || strcmp(options.occupations, 'insulator') % insulators
    options.occupations=-1;
  end
  if isscalar(options.occupations) && options.occupations>=0 % smearing
    calc=[ calc sprintf(', smearing=("gaussian",%g)', options.occupations/Ha) ]; % in Hartree
  end
  if options.nbands > 0
    % could not find it in NWChem
    % calc = [ calc sprintf(', nbands=%i', options.nbands) ];
  end
  if options.toldfe > 0
    calc = [ calc sprintf(', convergence={''energy'':%g}', options.toldfe) ];
  end
  if options.nsteps > 0
    calc = [ calc sprintf(', iterations=%i', options.nsteps) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  % to add: diagonalization
  
case {'QUANTUM','QE','ESPRESSO','QUANTUMESPRESSO','QUANTUM-ESPRESSO','PHON'}
  % best potentials for QE: SSSP http://materialscloud.org/sssp/
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  % ASE is installed. We use it to create a proper POSCAR file, then we call sqw_phon (QE)
  poscar = fullfile(options.target,'POSCAR_ASE');
  read = [ read 'from ase.io import write; ' ...
    'write("' poscar '",atoms, "vasp");' ...
    'write("' fullfile(target, 'configuration.png') '", atoms);' ...
    'write("' fullfile(target, 'configuration.eps') '", atoms);' ...
    'write("' fullfile(target, 'configuration.pov') '", atoms);' ...
    'write("' fullfile(target, 'configuration.cif') '", atoms, "cif");' ...
    'write("' fullfile(target, 'configuration.x3d') '", atoms, "x3d");' ...
    'write("' fullfile(target, 'configuration.pdb') '", atoms, "pdb");' ...
    'write("' fullfile(target, 'configuration.html') '", atoms, "html");' ...
    'write("' fullfile(target, 'configuration.etsf') '", atoms, "etsf");' ...
    'write("' fullfile(target, 'configuration_SHELX.res') '", atoms, "res");' ...
    'write("' fullfile(target, 'configuration_VASP') '", atoms, "vasp");' ];
  if isunix, setenv('LD_LIBRARY_PATH',''); end
  try
    [st, result] = system([ precmd 'python -c ''' read '''' ]);
  catch
    st = 127;
  end
  if isunix, setenv('LD_LIBRARY_PATH',ld_library_path); end
  if st ~= 0
    sqw_phonons_error([ mfilename ': failed converting input to POSCAR ' ...
      poscar ], options);
  end
% QE specific options:
%   options.mixing_ndim=scalar             number of iterations used in mixing
%     default=8. If you are tight with memory, you may reduce it to 4. A larger
%     value will improve the SCF convergence, and use more memory.
%   options.mixing_beta=scalar             mixing factor for self-consistency
%     default=0.7. use 0.3 to improve convergence
%   options.ecutrho=scalar                 kinetic energy cutoff (Ry) for charge
%     density and potential. Default=4*ecutwfc, suitable for PAW. Larger value
%     improves convergence, especially for ultra-soft PP (use 8-12*ecutwfc).
%   options.electron_maxstep=scalar        max number of iterations for SCF.
%     default=100. Larger value improves convergence.
%   options.conv_thr=scalar                Convergence threshold for 
%     selfconsistency. default=1e-6.
%   options.mpi    =scalar                 number of CPUs to use for PWSCF
%     this option requires MPI to be installed (e.g. openmpi).
%
  if options.nbands && ~isfield(options,'nbnd')
    options.nbnd=options.nbands; end
  if options.ecut > 0 && ~isfield(options,'ecutwfc')
    options.ecutwfc=options.ecut/Ry; end % in Ry
  if options.nsteps, options.electron_maxstep=options.nsteps; end
  if isscalar(options.occupations), options.occupations=options.occupations/Ry; end
  options.mpirun = status.mpirun;

  disp([ mfilename ': calling sqw_phon(' poscar ') with PHON/Quantum Espresso' ]);
  options.dos = 1;
  sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'init', options, [ 'sqw_phon(''' poscar ''', options); % QuantumEspresso/PHON wrapper' ]);
  signal=sqw_phon(poscar, options);


otherwise
  sqw_phonons_error([ mfilename ': Unknown calculator ' options.calculator ], options);
end
% ------------------------------------------------------------------------------
% end of specific parts for calculators














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
  if exist(configuration)
    copyfile(configuration, target);
  end
  
  % call python script with configuration
  read0 = [ read 'from ase.io import write; ' ...
    'write("' fullfile(target, 'configuration.png') '", atoms);' ...
    'write("' fullfile(target, 'configuration.eps') '", atoms);' ...
    'write("' fullfile(target, 'configuration.pov') '", atoms);' ...
    'write("' fullfile(target, 'configuration.cif') '", atoms, "cif");' ...
    'write("' fullfile(target, 'configuration.x3d') '", atoms, "x3d");' ...
    'write("' fullfile(target, 'configuration.pdb') '", atoms, "pdb");' ...
    'write("' fullfile(target, 'configuration.html') '", atoms, "html");' ...
    'write("' fullfile(target, 'configuration.etsf') '", atoms, "etsf");' ...
    'write("' fullfile(target, 'configuration_SHELX.res') '", atoms, "res");' ...
    'write("' fullfile(target, 'configuration_VASP') '", atoms, "vasp");' ];
  
  if isunix, setenv('LD_LIBRARY_PATH',''); end
  try
    [st, result] = system([ precmd 'python -c ''' read0 '''' ]);
  catch
    st = 127;
  end
  if isunix, setenv('LD_LIBRARY_PATH',ld_library_path); end
  if st ~= 0
    disp([ mfilename ': failed read input ' ...
      configuration ], options);
  end
  sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'init', options, calc);

  % call python script with calculator
  cd(target)
  disp([ mfilename ': creating Phonon/ASE model in ' target ]);
  disp([ '  ' configuration ]);
  disp([ '  ' calc ]);

  if options.gui && ishandle(options.gui), waitbar(0.15, options.gui, [ mfilename ': computing DynMatrix ASE/' options.calculator ' (be patient)' ]); end

  result = '';
  if isunix, setenv('LD_LIBRARY_PATH',''); end
  try
    [st, result] = system([ precmd 'python sqw_phonons_build.py' ]);
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
    sqw_phonons_error([ mfilename ': ' options.calculator ' failed. May be a convergence issue. Temporary files are in ' target ], options)
  end
  if exist(configuration)
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

  if exist(configuration)
    signal.UserData.configuration = fileread(configuration);
  else
    signal.UserData.configuration = configuration;
  end
  if isfield(options,'dos'), signal.UserData.DOS=[]; end
  signal.UserData.dir           = target;
  signal.UserData.options       = options;

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
    '  % import FREQ', ...
    '  FREQ=load(''FREQ'',''-ascii''); % in meV', ...
    'if isfield(this.UserData, ''DOS'') && isempty(this.UserData.DOS)', ...
    '  DOS = load(''DOS'',''-ascii''); DOS_w = load(''DOS_w'',''-ascii''); DOS=iData(DOS_w,DOS./sum(DOS));', ...
    '  DOS.Title = [ ''DOS '' strtok(this.Name) ]; xlabel(DOS,''DOS Energy [meV]'');', ...
    '  DOS.Error=0; this.UserData.DOS=DOS;', ...
    'end', ...
    'catch; disp([ ''model '' this.Name '' '' this.Tag '' could not run Python/ASE from '' target ]);', ...
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
disp([ '  in ' options.target ]);

if isfield(options, 'dos') && ~strcmpi(options.calculator, 'QUANTUMESPRESSO')
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
end
fprintf(1, '%s\n', cite{:});
disp('You can now evaluate the model using e.g.:')
disp('    qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);');
disp('    f=iData(s,[],qh,qk,ql,w); plot3(log(f(1,:, :,:)));');

sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'done', options, cite);

% handle autoplot option
if options.autoplot
  if options.gui && ishandle(options.gui), waitbar(0.75, options.gui, [ mfilename ': plotting phonons and DOS' ]); end
  f = sqw_phonons_plot(signal);
  if options.gui && ishandle(options.gui), delete(options.gui); end
else
  f = [];
end

sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'plot', options, f, signal);

% ------------------------------------------------------------------------------
function f = sqw_phonons_plot(signal)
  if isempty(signal) || ~isfield(signal.UserData, 'options') , return; end
  
  options = signal.UserData.options;
  disp([ mfilename ': Model ' options.configuration ' plotting phonons.' ])
  qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,150,151);
  fig=figure; 
  if options.dos, subplot(1,2,1); end
  f=iData(signal,[],qh,qk,ql,w);
  g=log(-f(1,:, :,:)); scatter3(g,'filled'); axis tight;
  view([38 26]);
  if options.dos
    subplot(1,2,2);
    try
    plot(signal.UserData.DOS); % plot the DOS, as indicated during model creation
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
  disp([ '<a href="matlab:doc(''' mfilename ''')">' mfilename ' help</a>' ])
end
sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'error', options, message);
error(message);

% ------------------------------------------------------------------------------



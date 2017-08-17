function [decl,calc,signal,options] = sqw_phonons_calc(options, status, calc_choice, read)
% sqw_phonons_calc: set python code to initiate the calculator
%   requires: nothing except options, status(requirements), choice and 
%   output:   python snippets to setup ASE Calculator, or signal=iFunc (QE case)

calc = ''; decl = ''; signal=[];
if nargin < 3, calc_choice = ''; end
if nargin < 4, read = ''; end

if isempty(calc_choice), calc_choice=options.calculator; end

% unit conversions
Ha = 27.2; Ry=Ha/2;

if isempty(status.(lower(options.calculator))) && isempty(options.command)
  sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  return
end

if     strncmp(options.occupations, 'auto',4) && ~strcmp(upper(calc_choice),'ELK')
  % Elk handle 'auto' by itself
  options.occupations=0.13;
elseif strncmp(options.occupations, 'semi',4)
  options.occupations=0.01;
elseif strncmp(options.occupations, 'metal',5)
  options.occupations=0.27;
elseif strncmp(options.occupations, 'fix',3) || strncmp(options.occupations, 'insu',4)
  options.occupations=-1;
end
  
switch upper(calc_choice)
% ==============================================================================
case 'ABINIT'
  % build: ./configure --enable-mpi CC=mpicc CXX=mpicxx FC=mpif90 --enable-optim --with-dft-flavor=libxc
  % ABINIT with pawxml requires ecut=1500 and nsteps > 30 (default)
  
  if isempty(strfind(status.(lower(options.calculator)),'abinis')) && isempty(options.command)
    options.command = status.(lower(options.calculator));
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ options.mpirun ' ' options.command ]; 
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
    elseif any(strcmpi(options.potentials, {'fhi', 'hgh', 'hgh.sc', 'hgh.k', 'tm', 'paw','pawxml'}))
      if strcmpi(options.potentials, 'paw') || strcmpi(options.potentials, 'pawxml')
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
  calc = 'calc = Abinit(chksymbreak=0, maxnsym=10384';
  if options.ecut <= 0, 
    options.ecut=1500; 
  end % no default in ABINIT (eV)
  if options.ecut > 0
    calc = [ calc sprintf(', ecut=%g', options.ecut/Ha) ];
    if isfield(options,'pps') && (strcmpi(options.pps, 'paw') || strcmpi(options.pps, 'pawxml'))
      if ~isfield(options,'pawecutdg')
        options.pawecutdg = 3*options.ecut;
      end
    end
  end
  if isfield(options,'pps') && strcmpi(options.pps, 'pawxml')
    options.ixc = 1; % the XC is stored in the PAWXML (libXC)
  end
  if isfield(options,'pawecutdg') && options.pawecutdg > 0
    calc = [ calc sprintf(', pawecutdg=%g', options.pawecutdg/Ha) ];
  end
  if isfield(options,'tolvrs') && options.tolvrs > 0
    calc = [ calc sprintf(', tolvrs=%g', options.tolvrs/Ha) ];
  else
    % if options.toldfe <= 0, options.toldfe=1e-6; end % in eV, necessary
    if options.toldfe > 0
      calc = [ calc sprintf(', toldfe=%g', options.toldfe/Ha) ];
    end
  end
  if ~isfield(options, 'pps') || isempty(options.pps)
    options.iscf=17;
  else
    calc = [ calc sprintf(', pps="%s"', options.pps) ];
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
    calc = [ calc sprintf(', xc="%s"', options.xc) ];
  end
  if isfield(options,'ixc')
    calc = [ calc sprintf(', ixc=%d', options.ixc) ];
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    % nbdblock, npband, AUTOPARAL=1
    % calc = [ calc sprintf(', nbdblock=%i', options.mpi) ];
    calc = [ calc sprintf(', autoparal=1') ];
  end
  if options.nbands > 0
    calc = [ calc sprintf(', nband=%i', options.nbands) ];
  end
  % the default nstep=30 is clearly not enough for ABINIT. Set to 100
  if options.nsteps <=100, options.nsteps=100; end
  if options.nsteps > 0
    calc = [ calc sprintf(', nstep=%i', options.nsteps) ];
  end
  
  if isscalar(options.occupations) && options.occupations >=0
    calc=[ calc sprintf(', tsmear=%g, occopt=4', options.occupations/Ha) ]; % "Cold smearing" of N. Marzari
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];

case 'ELK' % ===================================================================
  % requires custom compilation with elk/src/modmain.f90:289 maxsymcrys=1024

  % location of ELF pseudo-potentials is mandatory
  if isempty(options.potentials) && isempty(getenv('ELK_SPECIES_PATH'))
    if isunix, options.potentials = '/usr/share/elk-lapw/species';
      disp([ mfilename ': ' options.calculator ': assuming atom species are in' ])
      disp([ '  ' options.potentials ])
      disp('  WARNING: if this is not the right location, use options.potentials=<location>');
    else
      sqw_phonons_error([ mfilename ': ' options.calculator ': undefined "species". Use options.potentials=<location of elk/species>.' ], options)
      return
    end
  end
  if ~isempty(options.potentials)
    setenv('ELK_SPECIES_PATH', [ options.potentials, filesep ]);
  end
  if ~strcmp(status.(lower(options.calculator)),'elk') && isempty(options.command)
    options.command = status.(lower(options.calculator));
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ options.mpirun ' ' options.command ]; 
  end
  if ~isempty(options.command)
    cmd = options.command;
    if isempty(strfind(cmd, 'elk.out'))
      cmd = [ cmd ' > elk.out' ];
    end
    setenv('ASE_ELK_COMMAND', cmd);
  end
  
  decl = 'from ase.calculators.elk import ELK';
  calc = 'calc = ELK(tforce=True, tasks=0, mixtype=3'; % Pulay mixing

  if strcmp(options.occupations, 'auto')
    % stypen: 1-2: MethfesselPaxton; 3:FermiDirac
    calc=[ calc sprintf(', stype=3, autoswidth=True') ];
  elseif isscalar(options.occupations) && options.occupations >=0
    calc=[ calc sprintf(', stype=3, swidth=%g', options.occupations/Ha) ];
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=(%i,%i,%i)', options.kpoints) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc="%s"', options.xc) ];
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
  if isfield(options,' rgkmax') && options.rgkmax>0
    calc = [ calc sprintf(', rgkmax=%g', options.rgkmax) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  
case 'EMT'
  decl = 'from ase.calculators.emt import EMT';
  calc = 'calc  = EMT()';
case 'GPAW' % ==================================================================
  
  decl = 'from gpaw import GPAW, PW, FermiDirac';
  calc = [ 'calc = GPAW(symmetry={"point_group": False}, txt="' fullfile(options.target,'gpaw.log') '"' ]; % because small displacement breaks symmetry
  
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
    calc = [ calc sprintf(', mode="%s"', options.mode) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc="%s"', options.xc) ];
  end
  if ~isempty(options.diagonalization)
    if strncmpi(options.diagonalization, 'dav', 3) options.diagonalization='dav'; 
    elseif strcmpi(options.diagonalization, 'cg')  options.diagonalization='cg';
    else options.diagonalization='rmm-diis'; end
    calc = [ calc sprintf(', eigensolver="%s"', options.diagonalization) ];
  end
  if ~isempty(options.potentials)
    calc = [ calc sprintf(', setups="%s"', options.potentials) ];
  end
  if options.nbands > 0
    calc = [ calc sprintf(', nbands=%i', options.nbands) ];
  end
  if options.nsteps > 0
    calc = [ calc sprintf(', maxiter=%i', options.nsteps) ];
  end
  if options.toldfe > 0
    calc = [ calc sprintf(', convergence={"energy":%g}', options.toldfe) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  
case 'JACAPO' % ================================================================
  % fails with ASE 3.14
  % Requires to define variables under Ubuntu
  if isunix
    if isempty(options.potentials), options.potentials='/usr/share/dacapo-psp'; end
    if isempty(options.command) && (~isfield(options,'mpi') || isempty(options.mpi) || options.mpi <= 1)
      options.command   =status.(lower(options.calculator)); 
    end
  end
  if ~isempty(options.potentials)
    setenv('DACAPOPATH', options.potentials);
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.jacapo_mpi; end
    options.command = [ options.mpirun ' ' options.command ]; 
  end
  if ~isempty(options.command)
    setenv('DACAPOEXE_SERIAL', options.command); % DACAPOEXE_PARALLEL
  end
  
  decl = 'from ase.calculators.jacapo import Jacapo';
  calc = 'calc = Jacapo(symmetry=False';
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc="%s"', options.xc) ];
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
  if ~isempty(options.diagonalization)
    if strncmp(options.diagonalization, 'dav',3), options.diagonalization='eigsolve';
    elseif strcmpi(options.diagonalization, 'cg'), options.diagonalization='resmin';
    else options.diagonalization='rmm-diis'; end
    calc = [ calc sprintf(', electronic_minimization={"method":"%s"}', ...
             options.diagonalization) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];  
  % other options
  if options.toldfe > 0
    % this fails in ASE 3.14
    % calc = [ calc sprintf('; calc.SetConvergenceParameters(%g) ', options.toldfe) ];
  end
  

case 'NWCHEM' % ================================================================

  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ options.mpirun ' ' options.command ]; 
  end
  if ~isempty(options.command)
    cmd = options.command;
    if isempty(strfind(cmd, 'PREFIX.nw'))
      cmd = [ cmd ' PREFIX.nw > PREFIX.out' ];
    end
    setenv('ASE_NWCHEM_COMMAND', cmd);
  end
  
  
  
  decl = 'from ase.calculators.nwchem import NWChem';
  calc = sprintf('calc = NWChem(xc="%s", odft=True', ...
    options.xc);
  % check if we use KPTS
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', raw="""nwpw\\n  monkhorst-pack %i %i %i\\nend"""', options.kpoints) ];
  end
  if ~isempty(options.potentials)
    calc = [ calc sprintf(', basis="%s"', options.potentials) ];
  end
  % smearing is in Hartree
  if isscalar(options.occupations) && options.occupations>=0 % smearing
    calc=[ calc sprintf(', smearing=("gaussian",%g)', options.occupations/Ha) ]; % in Hartree
  end
  if (options.ecut > 0)
    % seems to be un-supported
    % calc = [ calc sprintf(', cutoff=%g', options.ecut) ];
  end
  if options.nbands > 0
    % could not find it in NWChem
    % calc = [ calc sprintf(', nbands=%i', options.nbands) ];
  end
  if options.toldfe > 0
    calc = [ calc sprintf(', convergence={"energy":%g}', options.toldfe) ];
  end
  if options.nsteps > 0
    calc = [ calc sprintf(', iterations=%i', options.nsteps) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  % to add: diagonalization
  
  
% ==============================================================================
case 'OCTOPUS'
  % experimental: Non-orthogonal cells support not implemented.
  
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ options.mpirun ' ' options.command ]; 
  end
  if isempty(options.command), options.command=status.(lower(options.calculator)); end
  if ~isempty(options.command)
    setenv('ASE_OCTOPUS_COMMAND', options.command);
  end
  
  decl = 'from ase.calculators.octopus import Octopus';
  calc = 'calc = Octopus(Output="dos + density + potential", OutputFormat="xcrysden", Spacing=0.25';
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', KPointsGrid=[[%i,%i,%i]], KPointsUseSymmetries=True', options.kpoints) ];
  end
  % pps: PseudopotentialSet can be:
  % standard: The standard set of Octopus that provides LDA pseudopotentials in the PSF format for some elements: H, Li, C, N, O, Na, Si, S, Ti, Se, Cd.
  % sg15: (experimental) The set of Optimized Norm-Conserving Vanderbilt PBE pseudopotentials (M. Schlipf and F. Gygi, Comp. Phys. Commun. doi:10.1016/j.cpc.2015.05.011 (2015)).
  % hgh_lda: The set of Hartwigsen-Goedecker-Hutter LDA pseudopotentials for elements from H to Rn. For many species a semi-core variant is available, it can be obtained by appending _sc to the element name. Ref: C. Hartwigsen, S. Goedecker, and J. Hutter, Phys. Rev. B 58, 3641 (1998).
  % hscv_lda: (experimental) The set of Hamann-Schlueter-Chiang-Vanderbilt (HSCV) potentials for LDA exchange and correlation downloaded from
  % hscv_pbe: (experimental) PBE version of the HSCV pseudopotentials. Check the documentation of the option hscv_lda for details and warnings.
  if ~isfield(options,'pps') || isempty(options.pps)
    if strcmpi(options.xc,'PBE') || strcmpi(options.potentials,'PBE')
      options.pps = 'pbe';
    elseif strcmpi(options.xc,'LDA') || strcmpi(options.potentials,'LDA')
      options.pps = 'lda';
    elseif strcmpi(options.potentials,'NC')
      options.pps = 'nc';
    elseif strcmpi(options.potentials,'hgh')
      options.pps = 'hgh';
    end
  end
  if isfield(options,'pps')
    if strcmpi(options.pps, 'pbe')
      options.pps = 'hscv_pbe';
    elseif strcmpi(options.pps, 'lda')
      options.pps = 'hscv_lda';
    elseif strcmpi(options.pps, 'hgh')
      options.pps = 'hgh_lda';
    elseif strcmpi(options.pps, 'nc')
      options.pps = 'sg15'; % ONCV
    end
    calc = [ calc sprintf(', ExperimentalFeatures=True, PseudopotentialSet="%s"', options.pps) ];
  end
  
  if isscalar(options.occupations) && options.occupations>=0 % smearing in eV, Gaussian
    calc=[ calc sprintf(', Smearing=%g, SmearingFunction="spline_smearing"', options.occupations) ];
  end
  
  if options.nsteps > 0
    calc = [ calc sprintf(', MaximumIter=%i', options.nsteps) ];
  end
  
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];


% ==============================================================================
case 'QUANTUMESPRESSO_ASE'
  
  
  decl = 'from qeutil import QuantumEspresso';
  calc = [ 'calc = QuantumEspresso(use_symmetry=False, tstress=True, nspin=2,  label=atoms.get_chemical_formula(), wdir="' options.target '"' ]; % because small displacement breaks symmetry
  
  if isscalar(options.occupations) && options.occupations>=0 % smearing in eV
    calc=[ calc sprintf(', occupations="smearing", smearing="methfessel-paxton", degauss=%g', options.occupations/Ry) ];
  elseif isscalar(options.occupations) && options.occupations < 0
    calc = [ calc sprintf(', occupations="fixed"') ]
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=(%i,%i,%i)', options.kpoints) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc="%s"', options.xc) ];
  end
  if ~isempty(options.potentials)
    calc = [ calc sprintf(', pseudo_dir="%s"', options.potentials) ];
  end
  if ~isempty(options.diagonalization)
    if strncmpi(options.diagonalization, 'dav', 3) options.diagonalization='david'; 
    else options.diagonalization='cg'; end
    calc = [ calc sprintf(', diagonalization="%s"', options.diagonalization) ];
    if ~isfield(options, 'mixing_beta'), options.mixing_beta = 0.7; end
  end
  if options.nbands > 0
    calc = [ calc sprintf(', nbnd=%i', options.nbands) ];
  end
  if isfield(options, 'toldfe') && ~isempty(options.toldfe)
    calc = [ calc sprintf(', conv_thr=%f', options.toldfe) ];
  end
  if isfield(options, 'mixing_beta') && ~isempty(options.mixing_beta)
    calc = [ calc sprintf(', mixing_beta=%f', options.mixing_beta) ];
  end
  if isfield(options, 'mixing_ndim') && ~isempty(options.mixing_ndim)
    calc = [ calc sprintf(', mixing_ndim=%i', options.mixing_ndim) ];
  end
  if isfield(options, 'mixing_beta') || isfield(options, 'mixing_ndim')
    calc = [ calc  ', mixing_mode = "plain"' ];
  end
  if options.nsteps>0
    calc = [ calc sprintf(', electron_maxstep=%i', options.nsteps) ];
  end
  if (options.ecut > 0)
    calc = [ calc sprintf(', ecutwfc=%g', options.ecut/Ry) ];
  else
    calc = [ calc ', ecutwfc=15*len(set(atoms.get_atomic_numbers()))' ];
  end
  if isfield(options, 'ecutrho') && ~isempty(options.ecutrho)
    calc = [ calc sprintf(', ecutrho=%f', options.ecutrho/Ry) ];
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    calc = [ calc sprintf(', procs=%i', options.mpi) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  

% ==============================================================================
case {'QUANTUMESPRESSO','QUANTUMESPRESSO_PHON'}
  % best potentials for QE: SSSP http://materialscloud.org/sssp/

  
  poscar = fullfile(options.target, 'configuration_POSCAR');  % this is a POSCAR file
  
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
%   options.toldfe=scalar                Convergence threshold for 
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

  decl = [ 'sqw_phon(''' poscar ''', options); % QuantumEspresso/PHON wrapper' ];
  % sqw_phonons_htmlreport(fullfile(options.target, 'sqw_phonons.html'), 'init', options, decl);
  
  signal=sqw_phon(poscar, options);
  if isempty(signal), decl=[]; return; end
  
  % get 'atoms' back from python
  if ~isempty(dir(fullfile(options.target, 'properties.mat')))
    signal.UserData.properties = load(fullfile(options.target, 'properties.mat'));
    if isfield(signal.UserData.properties, 'chemical_symbols')
      b_coh = sqw_phonons_b_coh(signal.UserData.properties.chemical_symbols);
      signal.UserData.properties.b_coh = b_coh;
    end
  else
    signal.UserData.properties = [];
  end
  signal.UserData.calc = 'quantumespresso';
  signal.UserData.configuration = fileread(poscar);
  
% ==============================================================================
case 'VASP'

  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ options.mpirun ' ' options.command ]; 
  end
  if isempty(options.command), options.command=status.(lower(options.calculator)); end
  if ~isempty(options.command)
    setenv('VASP_COMMAND', options.command);
  end
  if isunix
    if isempty(options.potentials), options.potentials='/usr/share/vasp/pseudo/'; end
  end
  if ~isempty(options.potentials) && isdir(options.potentials)
    setenv('VASP_PP_PATH', options.potentials);
  end

  decl = 'from ase.calculators.vasp import Vasp';
  calc = [ 'calc = Vasp(prec="Accurate", lreal=False, ibrion=-1, ' ...
           'nsw=0, lwave=False, lcharg=False, isym=2, ispin=2' ];
  
  % prec: Normal, Medium, Accurate -> predefined settings for encut
  % algo: Normal (Davidson) | Fast | Very_Fast (RMM-DIIS)
  % ibrion: -1: no ionic moves, -> nsw=0. 
  %   IBRION=8 computes full force constants in a single step.
  % ismear: -5 Blochl -4-tet -1-fermi 0-gaus >0 MP
  % setups: pv, sv 
  % nsw:  number of ionic steps. Default :	0
  % isym: switch symmetry stuff ON (1 or 2) or OFF (0)
  % could use: lepsilon = "T"
  if (options.ecut > 0)
    calc = [ calc sprintf(', encut=%g', options.ecut) ];
  end
  if options.toldfe <= 0, options.toldfe=1e-8; end % in eV, necessary
  if (options.toldfe > 0)
    calc = [ calc sprintf(', ediff=%g', options.toldfe) ];
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=[%i,%i,%i]', options.kpoints) ];
    if all(options.kpoints <= 1)
      calc = [ calc sprintf(', gamma=True') ];
    end
  end
  if ~isempty(options.diagonalization)
    if strncmp(options.diagonalization, 'dav',3), options.diagonalization='Normal';
    elseif strcmpi(options.diagonalization, 'cg'), options.diagonalization='Fast';
    else options.diagonalization='Very_Fast'; end
    calc = [ calc sprintf(', algo="%s"', options.diagonalization) ];
  end
    
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc="%s"', options.xc) ];
  end
  if isfield(options, 'pps') && ~isempty(options.pps)
    if ischar(options.pps)
      calc = [ calc sprintf(', setups="%s"', options.pps) ];
    elseif isstruct(options.pps)
      calc = [ calc sprintf(', setups={') ];
      f = fieldnames(options.pps);
      for index=1:numel(f)
        if ischar(options.pps.(f{index}))
          if index==1
            calc = [ calc sprintf('"%s":"%s"') ];
          else
            calc = [ calc sprintf(',"%s":"%s"') ];
          end
        end
      end
      calc = [ calc sprintf('}') ];
    end
  end
  if options.nbands > 0
    calc = [ calc sprintf(', nbands=%i', options.nbands) ];
  end
  if options.nsteps > 0
    calc = [ calc sprintf(', nelm=%i', options.nsteps) ];
  end
  
  if isscalar(options.occupations) && options.occupations >=0
    calc=[ calc sprintf(', sigma=%g, ismear=0', options.occupations) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  % this works FAST
  % calc=[ 'calc = Vasp(prec="medium", nsw=0, ediff=1e-6, nelm=60 ' ];
  
  calc = [ calc ')' ];


otherwise
  sqw_phonons_error([ mfilename ': Unknown calculator ' options.calculator ], options);
  return
end
% ------------------------------------------------------------------------------
% end of specific parts for calculators


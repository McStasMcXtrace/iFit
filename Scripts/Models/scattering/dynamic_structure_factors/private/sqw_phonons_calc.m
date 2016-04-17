function [decl,calc] = sqw_phonons_calc(options, status, calc_choice, read)

calc = ''; decl = '';
if nargin < 3, calc_choice = ''; end
if nargin < 4, read = ''; end

if isempty(calc_choice), calc_choice=options.calculator; end

Ha = 27.2; Ry=13.6;

switch upper(calc_choice)
% ==============================================================================
case 'ABINIT'
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  if isempty(strfind(status.(lower(options.calculator)),'abinis')) && isempty(options.command)
    options.command = status.(lower(options.calculator));
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
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
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
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
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
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
  calc = 'calc = GPAW(usesymm=False, txt="gpaw.log"'; % because small displacement breaks symmetry

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
    if strncmpi(options.diagonalization, 'dav', 3) options.diagonalization='dav'; 
    elseif strcmpi(options.diagonalization, 'cg')  options.diagonalization='cg';
    else options.diagonalization='rmm-diis'; end
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
    if isempty(options.command) && (~isfield(options,'mpi') || isempty(options.mpi) || options.mpi <= 1)
      options.command   =status.(lower(options.calculator)); 
    end
  end
  if ~isempty(options.potentials)
    setenv('DACAPOPATH', options.potentials);
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.jacapo_mpi; end
    options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
  end
  if ~isempty(options.command)
    setenv('DACAPOEXE_SERIAL', options.command); % DACAPOEXE_PARALLEL
  end
  
  decl = 'from ase.calculators.jacapo import Jacapo';
  calc = 'calc = Jacapo(symmetry=False';
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
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
    calc = [ calc sprintf('; calc.SetConvergenceParameters(%g) ', options.toldfe) ];
  end
  

case 'NWCHEM' % ================================================================
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
  end
  if ~isempty(options.command)
    cmd = options.command;
    if isempty(strfind(cmd, 'PREFIX.nw'))
      cmd = [ cmd ' PREFIX.nw > PREFIX.out' ];
    end
    setenv('ASE_NWCHEM_COMMAND', cmd);
  end
  
  
  
  decl = 'from ase.calculators.nwchem import NWChem';
  calc = sprintf('calc = NWChem(xc=''%s'', odft=True', ...
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
  if (options.ecut > 0)
    % seems to be un-supported
    % calc = [ calc sprintf(', cutoff=%g', options.ecut) ];
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
  
% ==============================================================================
case {'QUANTUM','QE','ESPRESSO','QUANTUMESPRESSO','QUANTUM-ESPRESSO','PHON'}
  % best potentials for QE: SSSP http://materialscloud.org/sssp/
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  
  poscar = fullfile(options.target,'POSCAR_ASE');
  if ismac,  precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else precmd=''; end
  try
    [st, result] = system([ precmd 'python -c ''' read '''' ]);
  catch
    st = 127;
  end
  if st ~= 0
      disp(read)
      disp(result)
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
  
  % get 'atoms' back from python
  if ~isempty(fullfile(target, 'atoms.mat'))
    signal.UserData.atoms = load(fullfile(target, 'atoms.mat'));
  else
    signal.UserData.atoms = [];
  end
  signal.UserData.calc = 'quantumespresso';
  signal.UserData.configuration = fileread(poscar);
  
% ==============================================================================
case 'VASP'
  if isempty(status.(lower(options.calculator))) && isempty(options.command)
    sqw_phonons_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  if isfield(options,'mpi') && ~isempty(options.mpi) && options.mpi > 1
    if isempty(options.command) options.command=status.(lower(options.calculator)); end
    options.command = [ status.mpirun ' -np ' num2str(options.mpi) ' ' options.command ]; 
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
  calc = 'calc = Vasp(isym=0, prec="Accurate", lreal="A", ibrion=2, nsw=5 ';
  % prec: Low, Normal, Accurate
  % algo: Normal (Davidson) | Fast | Very_Fast (RMM-DIIS)
  % ibrion
  % ismear: -5 Blochl -4-tet -1-fermi 0-gaus >0 MP
  % setups: pv, sv 
  if options.ecut <= 0, options.ecut=340; end % no default in ABINIT (eV)
  if (options.ecut > 0)
    calc = [ calc sprintf(', encut=%g', options.ecut) ];
  end
  if options.toldfe <= 0, options.toldfe=1e-5; end % in eV, necessary
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
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if isfield(options, 'pps') && ~isempty(options.pps)
    calc = [ calc sprintf(', setups=''%s''', options.pps) ];
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
  calc = [ calc ')' ];


otherwise
  sqw_phonons_error([ mfilename ': Unknown calculator ' options.calculator ], options);
end
% ------------------------------------------------------------------------------
% end of specific parts for calculators


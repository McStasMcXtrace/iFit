function signal=sqw_ph_ase(configuration, varargin)
% model=sqw_ph_ase(configuration, ..., options)
%
%   iFunc/sqw_ph_ase: computes phonon dispersions using the ASE.
%   A model which computes phonon dispersions from the forces acting between
%     atoms. The input argument is any configuration file describing the
%     material, e.g. CIF, PDB, POSCAR, ... supported by ASE.
%   The phonon spectra is computed using one of the calculator supported by the
%   Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase>.
%   Supported calculators are:
%     EMT       Effective Medium Theory calculator (Al,Cu,Ag,Au,Ni,Pd,Pt,H,C,N,O)
%     GPAW      Real-space/plane-wave/LCAO PAW code
%     NWChem    Gaussian based electronic structure code
%     Dacapo    Plane-wave ultra-soft pseudopotential code
%     ELK       Full Potential LAPW code
%     ABINIT    Plane-wave pseudopotential code
%
%   The simplest usage is to call: sqw_ph_ase('gui') which displays an entry dialog box
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
% options: an optional structure with optional settings:
%   options.target =path                   where to store all files and FORCES
%     a temporary directory is created when not specified.
%   options.supercell=scalar or [nx ny nz] supercell size. Default is 2.
%   options.calculator=string              EMT,GPAW,Elk,NWChem, Dacapo,ABINIT
%     Default=GPAW
%   options.dos=1                          options to compute the vibrational
%     density of states (vDOS) in UserData.DOS
%   options.potentials=string              basis datasets or pseudopotentials.
%     GPAW: see https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html#manual-setups
%     NWChem: see http://www.nwchem-sw.org/index.php/Release64:AvailableBasisSets
%     Dacapo, Elk, ABINIT: path to potentials.
%   options.kpoints=scalar or [nx ny nz]   Monkhorst-Pack grid
%   options.xc=string                      type of Exchange-Correlation functional to use
%     'LDA','PBE','revPBE','RPBE','PBE0','B3LYP'            for GPAW
%     'LDA','B3LYP','PBE','RHF','MP2'                       for NWChem
%     'LDA','PBE','REVPBE','PBESOL','WC06','AM05'           for ELK
%     'PZ','VWN','PW91','PBE','RPBE','revPBE'               for Dacapo/Jacapo
%     'LDA', 'PBE', 'revPBE', 'RPBE'                        for ABINIT
%     Default is 'PBE'.
%   options.mode='pw','fd', or 'lcao'      GPAW computation mode as Plane-Wave,
%     Finite Difference, or LCAO (linear combination of atomic orbitals). Default is 'fd'.
%   options.command='exe'                  Path to calculator executable
%   options.raw='name=value, ...'          Additional arguments passed to the ASE
%                                          using the Python syntax.
%
% options affecting memory usage:
%   options.diagonalization='dav' or 'cg' or 'rmm-diis' for GPAW
%     The Davidson method is faster than the conjugate-gradient but uses more
%     memory. The RMM-DIIS method allows parallelization.
%
% options affecting convergence:
%   options.occupations='metal'            for metals ('smearing') help converge
%                       'insulator'        for insulators
%                       'semiconductor'    sets 0 eV FermiDirac smearing
%                       'auto'             for Elk (automatic smearing)
%                       or 0 for semi-conductors
%                       or a value in eV for a FermiDirac distribution (GPAW,ABINIT)
%                                  in Hartree (NWChem,Elk)
%
%   options.ecutwfc=scalar                 kinetic energy cutoff (eV) for
%     wavefunctions. Default is 340 eV. Larger value improves convergence.
%
% The options can also be entered as a strings with 'field=value; ...'.
%
% Once the model has been created, its use requires that axes are given on
% regular qx,qy,qz grids.
%     
% Example:
%   s=sqw_ph_ase('bulk("Cu", "fcc", a=3.6, cubic=True)','EMT','metal');
%   qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,50,51);
%   f=iData(s,[],qh,qk,ql,w); scatter3(log(f(1,:, :,:)),'filled');
%   figure; plot(s.UserData.DOS); % plot the DOS, as indicated during model creation
%
%   s=sqw_ph_ase('bulk("Si", "diamond", a=5.4, cubic=True)','semiconductor');
%
%   s=sqw_ph_ase([ ifitpath 'Data/POSCAR_Al'],'dos','metal','EMT');
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
%
% input:  p: sqw_ph_ase model parameters (double)
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
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phon, sqw_cubic_monoatomic, sqw_sine3d, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>

persistent status

signal = [];

options= sqw_ph_ase_argin(varargin{:});

if nargin == 1 && strcmp(configuration,'gui')
  options.gui=1;
  configuration=[];
end

if nargin == 0 || isempty(configuration)
  configuration = 'bulk("Al", "fcc", a=4.05)';
end

if ~exist('status') || isempty(status) || ~isstruct(status)
  status = sqw_ph_ase_requirements;
end

if options.gui
  % pop-up a simple dialog requesting for:
  %  * configuration
  %  * calculator
  %  * metal, insulator, semiconductor
  %  * supercell
  %  * kpoints
  % and will select autoplot, 'dos'
  calcs = 'EMT';
  for index={'gpaw','elk','jacapo','nwchem','abinit'};
    if status.(index{1}), calcs = [ calcs ', ' upper(index{1}) ]; end
  end
  calcs = [ calcs ', QuantumEspresso' ];
  NL = sprintf('\n');
  prompt = { [ '{\bf Atom/molecule/system configuration}' NL 'a CIF/PDB/POSCAR/... name or e.g. bulk("Cu", "fcc", a=3.6, cubic=True),' NL  'molecule("H2O"), or nanotube(6, 0, length=4)' ], ...
  [ '{\bf Calculator}' NL 'one of ' calcs ], ...
  [ '{\bf Smearing}' NL 'metal, semiconductor, insulator or left empty' ], ...
  [ '{\bf Supercell}' NL 'the size of the repeated model = system*supercell (3-vector)' ], ...
  [ '{\bf K-Points}' NL 'Monkhorst-Pack grid which determines the K sampling (3-vector)' ] };
  dlg_title = 'iFit: Model: phonons from ASE';
  defAns    = {configuration, options.calculator, options.occupations, ...
    num2str(options.supercell), num2str(options.kpoints)};
  num_lines = [ 1 1 1 1 1 ]';
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
  options.supercell  = str2num(answer{4});
  options.kpoints    = str2num(answer{5});
  options.autoplot   = 1;
  options.dos        = 1;
  options.gui        = waitbar(0, [ mfilename ': starting' ], 'Name', [ 'iFit: ' mfilename ' ' configuration ]);
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

% handle supported calculators: decl and calc
if options.gui && ishandle(options.gui), waitbar(0.05, options.gui, [ mfilename ': configuring ' options.calculator ]); end
if isunix, precmd = 'LD_LIBRARY_PATH= ; '; else precmd=''; end
switch upper(options.calculator)
case {'QUANTUM','QE','ESPRESSO','QUANTUMESPRESSO','QUANTUM-ESPRESSO','PHON'}
  % ASE is installed. We use it to create a proper POSCAR file, then we call sqw_phon (QE)
  poscar = fullfile(options.target,'POSCAR_ASE');
  read = [ read 'from ase.io import write; write("' poscar '",atoms, format="vasp")' ];
  try
    [st, result] = system([ precmd 'python -c ''' read '''' ]);
  catch
    st = 127;
  end
  if st ~= 0
    sqw_ph_ase_error([ mfilename ': failed calling sqw_phon with ' ...
      poscar ], options);
  end
  if options.gui && ishandle(options.gui), 
    delete(options.gui);
  end
  options.gui = 0;
  disp([ mfilename ': calling sqw_phon(' poscar ') with PHON/Quantum Espresso' ]);
  signal=sqw_phon(poscar, options);
  return
case 'ABINIT'
  if ~status.(lower(options.calculator)) && isempty(options.command)
    sqw_ph_ase_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  if ~isempty(options.command)
    cmd = options.command;
    if isempty(strfind(cmd, 'PREFIX.files'))
      cmd = [ cmd ' < PREFIX.files > PREFIX.log' ];
    end
    setenv('ASE_ABINIT_COMMAND', cmd);
  end
  if ~isempty(options.potentials)
    setenv('ABINIT_PP_PATH', options.potentials)
  end
  decl = 'from ase.calculators.abinit import Abinit';
  calc = 'calc = Abinit(chksymbreak=0, toldfe=1.0e-5';
  if options.ecutwfc <= 0, options.ecutwfc=200; end % no default in ABINIT
  if (options.ecutwfc > 0)
    calc = [ calc sprintf(', ecut=%g', options.ecutwfc) ];
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=[%i,%i,%i]', options.kpoints) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
case 'ELK' % ===================================================================
  % requires custom compilation with elk/src/modmain.f90:289 maxsymcrys=1024
  if ~status.(lower(options.calculator)) && isempty(options.command)
    sqw_ph_ase_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  % location of ELF pseudo-potentials is mandatory
  if isempty(options.potentials) && isempty(getenv('ELK_SPECIES_PATH'))
    if isunix, options.potentials = '/usr/share/elk-lapw/species';
      disp([ mfilename ': ' options.calculator ': assuming atom species are in' ])
      disp([ '  ' options.potentials ])
      disp('  WARNING: if this is not the right location, use options.potentials=<location>');
    else
      sqw_ph_ase_error([ mfilename ': ' options.calculator ': undefined "species". Use options.potentials=<location of elk/species>.' ], options)
    end
  end
  if ~isempty(options.potentials)
    setenv('ELK_SPECIES_PATH', [ options.potentials, filesep ]);
  end
  % test if ELK executable is named elk-lapw
  [st,result]=system([ precmd 'elk-lapw' ]);
  if isempty(options.command) && st == 0
    options.command = 'elk-lapw';
  end
  if ~isempty(options.command)
    cmd = options.command;
    if isempty(strfind(cmd, 'elk.out'))
      cmd = [ cmd ' > elk.out' ];
    end
    setenv('ASE_ELK_COMMAND', cmd);
  end
  
  decl = 'from ase.calculators.elk import ELK';
  calc = 'calc = ELK(tforce=True, tasks=0, rgkmax=4.0, mixtype=3';
  if strcmp(options.occupations, 'smearing') || strcmp(options.occupations, 'metal') % metals
    options.occupations=0.1;
  elseif strcmp(options.occupations, 'semiconductor')
    options.occupations=0;
  elseif strcmp(options.occupations, 'fixed') || strcmp(options.occupations, 'insulator') % insulators
    options.occupations=-1;
  end
  if strcmp(options.occupations, 'auto')
    % other distribution: MethfesselPaxton
    calc=[ calc sprintf(', stype=3, autoswidth=True') ];
  elseif isscalar(options.occupations) && options.occupations >=0
    calc=[ calc sprintf(', stype=3, swidth=%g', options.occupations) ];
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=(%i,%i,%i)', options.kpoints) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  
case 'EMT'
  decl = 'from ase.calculators.emt import EMT';
  calc = 'calc  = EMT()';
case 'GPAW' % ==================================================================
  if ~status.(lower(options.calculator)) && isempty(options.command)
    sqw_ph_ase_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  
  decl = 'from gpaw import GPAW, PW, FermiDirac';
  calc = 'calc = GPAW(usesymm=False'; % because small displacement breaks symmetry
  if strcmp(options.occupations, 'smearing') || strcmp(options.occupations, 'metal') % metals
    options.occupations=0.1;
  elseif strcmp(options.occupations, 'semiconductor')
    options.occupations=0;
  elseif strcmp(options.occupations, 'fixed') || strcmp(options.occupations, 'insulator') % insulators
    options.occupations=-1;
  end
  if isscalar(options.occupations) && options.occupations>=0 % smearing
    calc=[ calc sprintf(', occupations=FermiDirac(%g)', options.occupations) ];
    % other distribution: MethfesselPaxton
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=(%i,%i,%i)', options.kpoints) ];
  end
  if (options.ecutwfc > 0)
    options.mode=sprintf('PW(%g)', options.ecutwfc);
  end
  if ~isempty(options.mode)
    calc = [ calc sprintf(', mode=''%s''', options.mode) ];
  end
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if ~isempty(options.diagonalization)
    calc = [ calc sprintf(', eigensolver=''%s''', options.diagonalization) ];
  end
  if ~isempty(options.potentials)
    calc = [ calc sprintf(', setups=''%s''', options.potentials) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];
  
case 'JACAPO' % ================================================================
  if ~status.(lower(options.calculator)) && isempty(options.command)
    sqw_ph_ase_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
  end
  
  % Requires to define variables under Ubuntu
  if isunix
    if isempty(options.potentials), options.potentials='/usr/share/dacapo-psp'; end
    if isempty(options.command),    options.command   ='dacapo_serial.run'; end
  end
  if ~isempty(options.potentials)
    setenv('DACAPOPATH', options.potentials);
  end
  if ~isempty(options.command)
    setenv('DACAPOEXE_SERIAL', options.command);
  end
  
  decl = 'from ase.calculators.jacapo import Jacapo';
  calc = 'calc = Jacapo(symmetry=False';
  if ~isempty(options.xc)
    calc = [ calc sprintf(', xc=''%s''', options.xc) ];
  end
  if (options.ecutwfc <= 0)
    options.ecutwfc = 340;
  end
  if all(options.kpoints > 0)
    calc = [ calc sprintf(', kpts=(%i,%i,%i)', options.kpoints) ];
  end
  if (options.ecutwfc > 0)
    calc = [ calc sprintf(', pw=%g', options.ecutwfc) ];
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];  
  
case 'NWCHEM' % ================================================================
  if ~status.(lower(options.calculator)) && isempty(options.command)
    sqw_ph_ase_error([ mfilename ': ' options.calculator ' not available. Check installation' ], options)
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
    calc=[ calc sprintf(', smearing=("gaussian",%g)', options.occupations) ]; % in Hartree
  end
  if ~isempty(options.raw)
    calc = [ calc sprintf(', %s', options.raw) ];
  end
  calc = [ calc ')' ];

otherwise
  sqw_ph_ase_error([ mfilename ': Unknown ASE calculator ' options.calculator ], options);
end

if options.gui && ishandle(options.gui), waitbar(0.10, options.gui, [ mfilename ': writing python script ASE/' options.calculator ]); end

if strcmp(upper(options.calculator), 'GPAW')
  % GPAW Bug: gpaw.aseinterface.GPAW does not support pickle export for 'input_parameters'
  sav = sprintf('ph.calc=None\npickle.dump(ph, fid)');
else
  sav = 'pickle.dump(ph, fid)';
end
% start python --------------------------
script = { ...
  decl, ...
  'from ase.phonons import Phonons', ...
  'import pickle', ...
  '# Setup crystal and calculator', ...
  read, ...
  'from ase.io import write', ...
  'write("configuration.png", atoms)', ...
  'write("configuration.eps", atoms)', ...
  'write("configuration.pov", atoms)', ...
  calc, ...
  '# Phonon calculator', ...
sprintf('ph = Phonons(atoms, calc, supercell=(%i, %i, %i), delta=0.05)',options.supercell), ...
  'ph.run()', ...
  '# Read forces and assemble the dynamical matrix', ...
  'ph.read(acoustic=True)', ...
  '# save ph', ...
[ 'fid = open(''' target '/ph.pkl'',''wb'')' ], ...
  sav, ...
  'fid.close()' };
% end   python --------------------------

% write the script in the target directory
fid = fopen(fullfile(target,'sqw_ph_ase_build.py'),'w');
fprintf(fid, '%s\n', script{:});
fclose(fid);
% copy the configuration into the target
if exist(configuration)
  copyfile(configuration, target);
end

% call python script
cd(target)
disp([ mfilename ': creating Phonon/ASE model from ' target ]);
disp([ '  ' configuration ]);
disp([ '  ' calc ]);

if options.gui && ishandle(options.gui), waitbar(0.15, options.gui, [ mfilename ': computing DynMatrix ASE/' options.calculator ' (be patient)' ]); end

result = '';
try
  [st, result] = system([ precmd 'python sqw_ph_ase_build.py' ]);
  disp(result)
catch
  disp(result)
  sqw_ph_ase_error([ mfilename ': failed calling ASE with script ' ...
    fullfile(target,'sqw_ph_ase_build.py') ], options);
end
cd(pw)

% then read the pickle file to store it into the model
if options.gui && ishandle(options.gui), waitbar(0.7, options.gui, [ mfilename ': building model' ]); end
try
  signal.UserData.ph_ase = fileread(fullfile(target, 'ph.pkl')); % binary
catch
  
  sqw_ph_ase_error([ mfilename ': ' options.calculator ' failed. May be a convergence issue. Temporary files are in ' target ], options)
end
if exist(configuration)
  [dummy, signal.UserData.input]= fileparts(configuration);
else
  signal.UserData.input = configuration;
end

signal.Name           = [ 'Sqw_ASE_' signal.UserData.input ' Phonon/ASE DHO [' mfilename ']' ];

signal.Description    = [ 'S(q,w) 3D dispersion Phonon/ASE with DHO line shape. ' configuration ];

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
[ '  fid=fopen(fullfile(target,''sqw_ph_ase_eval.py''),''w'');' ], ...
[ '  fprintf(fid, ''# Script file for Python/ASE to compute the modes from ' configuration ' in %s\n'', target);' ], ...
  '  fprintf(fid, ''#   ASE: S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002\n'');', ...
  '  fprintf(fid, ''#   <https://wiki.fysik.dtu.dk/ase>\n'');', ...
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
[ '  [status,result] = system(''' precmd 'python sqw_ph_ase_eval.py'');' ], ...
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

% when model is successfully built, display citations for ASE
disp([ mfilename ': Model ' configuration ' built using: (please cite)' ])
disp([ '  in ' options.target ]);
if isfield(options, 'dos')
  disp('INFO: The vibrational density of states (vDOS) will be computed at first model evaluation.');
end
disp(' * Atomic Simulation Environment')
disp('     S. R. Bahn and K. W. Jacobsen, Comput. Sci. Eng., Vol. 4, 56-66, 2002.')
disp('     <https://wiki.fysik.dtu.dk/ase>. ')
disp(' * iFit:   E. Farhi et al, J. Neut. Res., 17 (2013) 5.')
disp('     <http://ifit.mccode.org>.')
switch upper(options.calculator)
case 'GPAW'
disp(' * GPAW:   J. J. Mortensen et al, Phys Rev B, Vol. 71, 035109 (2005).');
case 'NWCHEM'
disp(' * NWChem: M. Valiev et al, Comput. Phys. Commun. 181, 1477 (2010).');
case 'ELK'
disp(' * ELK:    http://elk.sourceforge.net');
case {'DACAPO','JACAPO'}
disp(' * DACAPO: B. Hammer et al, Phys. Rev. B 59, 7413 (1999).');
case 'ABINIT'
disp(' * ABINIT: X. Gonze et al, Computer Physics Communications 180, 2582-2615 (2009).');
case 'EMT'
disp(' * EMT:    K.W. Jacobsen et al, Surf. Sci. 366, 394–402 (1996).');
end

% handle autoplot option
if options.autoplot
  if options.gui && ishandle(options.gui), waitbar(0.75, options.gui, [ mfilename ': plotting phonons and DOS' ]); end
  disp([ mfilename ': Model ' configuration ' plotting phonons.' ])
  qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,150,151);
  figure; 
  if options.dos, subplot(1,2,1); end
  f=iData(signal,[],qh,qk,ql,w); scatter3(log(f(1,:, :,:)),'filled'); axis tight;
  view([38 26]);
  if options.dos
    subplot(1,2,2);
    plot(signal.UserData.DOS); % plot the DOS, as indicated during model creation
  end
  drawnow
  if options.gui && ishandle(options.gui), delete(options.gui); end
end










% ------------------------------------------------------------------------------
function status = sqw_ph_ase_requirements
% sqw_ph_ase_requirements: check for availability of ASE and MD codes
%
% returns a structure with a field for each MD software being 1 when available.
status = [];

% test for ASE in Python
if isunix, precmd = 'LD_LIBRARY_PATH= ; '; else precmd=''; end
[status.ase, result] = system([ precmd 'python -c "import ase.version; print ase.version.version"' ]);
if status.ase ~= 0
  disp([ mfilename ': error: requires ASE to be installed.' ])
  disp('  Get it at <https://wiki.fysik.dtu.dk/ase>.');
  disp('  Packages exist for Debian/Mint/Ubuntu, RedHat/Fedora/SuSE, MacOSX and Windows.');
  error([ mfilename ': ASE not installed' ]);
else
  disp([ mfilename ': using ASE ' result ]);
  disp('Available calculators:');
  status.emt=1;
  disp('  EMT           only for Al,Cu,Ag,Au,Ni,Pd,Pt,H,C,N,O');
  % test for GPAW
  [st, result] = system([ precmd 'python -c "from gpaw import GPAW"' ]);
  if st == 0
    status.gpaw=1;
    disp('  GPAW (http://wiki.fysik.dtu.dk/gpaw)');
  else
    status.gpaw=0;
  end
  % test for NWChem
  [st, result] = system([ precmd 'python -c "from ase.calculators.nwchem import NWChem"' ]);
  status.nwchem=0;
  if st == 0
    % now test executable
    % create a fake nwchem.nw
    f = tempname;
    dlmwrite([ f '.nw' ], '');
    [st,result]=system([ precmd 'nwchem ' f '.nw' ]);
    if st==0 || st==139
      status.nwchem=1;
      disp('  NWChem (http://www.nwchem-sw.org/) as "nwchem"');
    end
    delete([ f '.*' ])
    [p,f] = fileparts(f);
    delete([ f '.db' ])
  end
  
  % test for Jacapo
  [st, result] = system([ precmd 'python -c "from ase.calculators.jacapo import Jacapo"' ]);
  status.jacapo=0;
  if st == 0
    % now test executable
    [st,result]=system([ precmd 'dacapo_serial.run' ]);
    if st == 0
      status.jacapo=1;
      disp('  Dacapo (http://wiki.fysik.dtu.dk/dacapo) as "dacapo_serial.run"');
    end
  end
  % test for Elk
  [st, result] = system([ precmd 'python -c "from ase.calculators.elk import ELK"' ]);
  status.elk=0;
  if st == 0
    % now test executable
    [st,result]=system([ precmd 'elk' ]);
    if st ~= 0
      % try elk-lapw
      [st,result]=system([ precmd 'elk-lapw' ]);
      if st == 0
        status.elk=1;
        disp('  Elk (http://elk.sourceforge.net) as "elk-lapw"');
      end
    else
      status.elk=1;
      disp('  Elk (http://elk.sourceforge.net) as "elk"');
    end
    % test for ABINIT
    [st, result] = system([ precmd 'python -c "from ase.calculators.abinit import Abinit"' ]);
    status.abinit=0;
    if st == 0
      % now test executable
      [st,result]=system([ precmd 'echo "0" | abinis' ]);
      if st == 0 || st == 2
        status.abinit=1;
        disp('  ABINIT (http://www.abinit.org/) as "abinis"');
      end
    end
  end
  disp('Calculator executables can be specified as ''options.command=exe'' when building a model.');
  
  %  lj (lenard-jones)
  %  morse
  %  eam

end

% ------------------------------------------------------------------------------
function options=sqw_ph_ase_argin(varargin)
% sqw_ph_ase_argin: extracts options from the arguments
%
% returns an 'options' structure.

% defaults
options.supercell  = 2;
options.calculator = 'GPAW';
options.kpoints    = 1;
options.xc         = 'PBE';
options.mode       = 'fd';            % GPAW
options.potentials = '';
options.diagonalization = ''; % GPAW would prefer rmm-diis
options.occupations= '';
options.ecutwfc    = 0;
options.command    = '';
options.raw        = '';
options.autoplot   = 0;
options.gui        = 0;

% read input arguments
  for index=1:numel(varargin)
  if ischar(varargin{index}) && isempty(dir(varargin{index})) && ~isempty(find(varargin{index} == '='))
    % first try to build a structure from the string
    this = str2struct(varargin{index});
    if isstruct(this)
      varargin{index} = this;
    end
  end
  if ischar(varargin{index})
    [p,f,e] = fileparts(varargin{index});
    % handle static options: metal,insulator, random
    if strcmp(varargin{index},'smearing') || strcmp(varargin{index},'metal')
      options.occupations = 'smearing';
    elseif strcmp(varargin{index},'fixed') || strcmp(varargin{index},'insulator')
      options.occupations = 'fixed';
    elseif strcmp(varargin{index},'semiconductor')
      options.occupations = 'semiconductor';
    elseif strcmp(lower(varargin{index}),'dos')
      options.dos = 1;
    elseif strcmp(lower(varargin{index}),'emt')
      options.calculator = 'EMT';
    elseif strcmp(lower(varargin{index}),'gpaw')
      options.calculator = 'GPAW';
    elseif strcmp(lower(varargin{index}),'jacapo') || strcmp(lower(varargin{index}),'dacapo')
      options.calculator = 'Jacapo';
    elseif strcmp(lower(varargin{index}),'nwchem')
      options.calculator = 'NWChem';
    elseif strcmp(lower(varargin{index}),'elk')
      options.calculator = 'Elk';
    elseif strcmp(lower(varargin{index}),'abinit')
      options.calculator = 'ABINIT';
    elseif strcmp(lower(varargin{index}),'autoplot')
      options.autoplot = 1;
    elseif strcmp(lower(varargin{index}),'gui')
      options.gui = 1;
    end
  elseif isstruct(varargin{index})
    % a structure: we copy the fields into options.
    this = varargin{index};
    f    =fieldnames(this);
    for i=1:numel(fieldnames(this))
      options.(f{i}) = this.(f{i});
    end
  end
end
if ~isfield(options,'target')
  options.target = tempname; % everything will go there
  mkdir(options.target)
end

if isscalar(options.supercell)
  options.supercell=[ options.supercell options.supercell options.supercell ]; 
end
if isscalar(options.kpoints)
  options.kpoints=[ options.kpoints options.kpoints options.kpoints ]; 
end

% ------------------------------------------------------------------------------

function sqw_ph_ase_error(message, options)

if options.gui && ishandle(options.gui)
  delete(options.gui);
  errordlg(message, [ 'iFit: ' mfilename ' ' options.configuration ' FAILED' ]);
end
error(message);

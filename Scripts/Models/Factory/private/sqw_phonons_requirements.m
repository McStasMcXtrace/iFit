function status = sqw_phonons_requirements
% sqw_phonons_requirements: check for availability of ASE and MD codes
%
% MPI, EMT, GPAW, NWChem, Dacapo, Abinit, Elk, QE, VASP
%
% returns a structure with a field for each MD software being 1 when available.
status = [];

if ismac,  precmd = 'DYLD_LIBRARY_PATH= ;';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
else precmd=''; end

disp('Available packages:');

% test for python
status.python = '';
for calc={'python'}
  % now test executable
  [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
  if any(st == 0:2)
      status.python=calc{1};
      st = 0;
      disp([ '  Python          (http://www.python.org) as "' status.python '"' ]);
      break;
  end
end
if isempty(status.python)
  error([ mfilename ': Python not installed. This is required.' ]);
end

% test for ASE in Python

[status.ase, result] = system([ precmd status.python ' -c "import ase"' ]);
if status.ase ~= 0
  disp([ mfilename ': error: requires ASE to be installed.' ])
  disp('  Get it at <https://wiki.fysik.dtu.dk/ase>.');
  disp('  Packages exist for Debian/Mint/Ubuntu, RedHat/Fedora/SuSE, MacOSX and Windows.');
  error([ mfilename ': ASE not installed. This is required.' ]);
else
  disp([ mfilename ': using ASE ' result ]);
  status.emt='ase-run';
  status.ase=sscanf(result,'%d.%d');
  disp('  EMT             only for Al,Cu,Ag,Au,Ni,Pd,Pt,H,C,N,O');
  
  % test for mpirun
  status.mpirun = '';
  for calc={'mpirun','mpiexec'}
    % now test executable
    [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
    if any(st == 0:2)
        status.mpirun=calc{1};
        st = 0;
        disp([ '  MPI             (http://www.openmpi.org) as "' status.mpirun '"' ]);
        break;
    end
  end
  
  % test for PhonoPy
  [status.phonopy, result] = system([ precmd status.python ' -c "from phonopy import Phonopy"' ]);
  if any(st == 0)
    status.phonopy = 'phonopy';
    disp([ '  PhonoPy         (https://atztogo.github.io) as "' status.phonopy '"' ]);
  else
    status.phonopy = '';
  end
  
  % test for GPAW
  [st, result] = system([ precmd status.python ' -c "from gpaw import GPAW"' ]);
  if any(st == 0:2)
    status.gpaw='gpaw-python';
    disp([ '  GPAW            (http://wiki.fysik.dtu.dk/gpaw) as "' status.gpaw '"' ]);
  else
    status.gpaw='';
  end
  
  % test for NWChem
  [st, result] = system([ precmd status.python ' -c "from ase.calculators.nwchem import NWChem"' ]);
  status.nwchem='';
  if any(st == 0:2)
    % now test executable
    % create a fake nwchem.nw
    f = tempname;
    dlmwrite([ f '.nw' ], '');
    [st,result]=system([ precmd 'nwchem ' f '.nw' ]);
    if any(st == 0:2) || st==139
      status.nwchem='nwchem';
    end
    delete([ f '.*' ])
    [p,f] = fileparts(f);
    if ~isempty(dir([ f '.db' ])), delete([ f '.db' ]); end
  end
  if ~isempty(status.nwchem)
    disp(['  NWChem          (http://www.nwchem-sw.org/) as "' status.nwchem '"' ]);
  end
  
  % test for Jacapo
  [st, result] = system([ precmd status.python ' -c "from ase.calculators.jacapo import Jacapo"' ]);
  status.jacapo='';
  if st == 0
    % now test executable: serial
    for calc={'dacapo_serial.run','dacapo.run','dacapo'}
      [st,result]=system([ precmd calc{1} ]);
      if any(st == 0:2)
          status.jacapo=calc{1};
          st = 0;
          break;
      end
    end
    % now test executable: mpi
    [st,result]=system([ precmd 'dacapo_mpi.run' ]);
    if st == 0 || st == 2
      status.jacapo_mpi='dacapo_mpi.run';
    else
      status.jacapo_mpi='';
    end
  end
  if ~isempty(status.jacapo)
    if ~isempty(status.jacapo_mpi)
      disp([ '  Dacapo          (http://wiki.fysik.dtu.dk/dacapo) as "' status.jacapo '" and "' status.jacapo_mpi '"' ]);
    else
      disp([ '  Dacapo          (http://wiki.fysik.dtu.dk/dacapo) as "' status.jacapo '"' ]);
    end
  end
  
  % test for Elk
  [st, result] = system([ precmd status.python ' -c "from ase.calculators.elk import ELK"' ]);
  status.elk='';
  if any(st == 0:2)
    % now test executable
    for calc={'elk','elk-lapw'}
      [st,result]=system([ precmd calc{1} ]);
      if st == 0 || st == 2
          status.elk=calc{1};
          st = 0;
          break;
      end
    end
  end
  if ~isempty(status.elk)
    disp([ '  Elk             (http://elk.sourceforge.net) as "' status.elk '"' ]);
  end
  
  % test for ABINIT
  [st, result] = system([ precmd status.python ' -c "from ase.calculators.abinit import Abinit"' ]);
  status.abinit='';
  if any(st == 0:2)
    for calc={'abinit','abinis','abinip'}
      % now test executable
      [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
      if any(st == 0:2)
          status.abinit=calc{1};
          st = 0;
          break;
      end
    end
  end
  if ~isempty(status.abinit)
    disp([ '  ABINIT          (http://www.abinit.org/) as "' status.abinit '"' ]);
  end
  
  % test for QuantumEspresso
  status.quantumespresso = '';
  for calc={'pw.x','pw.exe','pw','pwscf'}
    % now test executable
    [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
    if any(st == 0:2)
        status.quantumespresso=calc{1};
        st = 0;
        break;
    end
  end
  
  % test for QE/ASE
  status.quantumespresso_ase = '';
  if ~isempty(status.quantumespresso)
    [st, result] = system([ precmd status.python ' -c "from qeutil import QuantumEspresso"' ]);
    if any(st == 0:2)
      status.quantumespresso_ase='qeutil';
      disp([ '  QEutil          (https://jochym.github.io/qe-doc/) as "' status.quantumespresso_ase '"' ]);
    else
      status.quantumespresso_ase='';
    end
  end
  
  % test for PHON
  [st, result] = system([ precmd 'phon' ]);
  try
    delete('CRASH');
    delete('input_tmp.in');
  end
  if any(st == 0:2)
    status.phon = 'phon';
  else
    status.phon = '';
    if isempty(status.quantumespresso_ase)
      status.quantumespresso = '';  % no PHON, nor QEutil
    end
  end
  if ~isempty(status.quantumespresso)
    disp([ '  QuantumEspresso (http://www.quantum-espresso.org/) as "' status.quantumespresso '"' ]);
  end
  
  % test for VASP
  [st, result] = system([ precmd 'vasp' ]);
  if any(st == 0:2)
    status.vasp = 'vasp';
  else
    status.vasp = '';
  end
  if ~isempty(status.vasp)
    disp([ '  VASP            (http://www.vasp.at/) as "' status.vasp '"' ]);
  end
  % must set:
  % VASP_COMMAND=vasp
  % VASP_PP_PATH=/usr/share/vasp/pseudo
  
  % test for Octopus
  status.octopus = '';
  for calc={'octopus','octopus_mpi'}
    % now test executable
    [st,result]=system([ precmd calc{1} ' -v' ]);
    if any(st == 0:2)
        status.octopus=calc{1};
        st = 0;
        break;
    end
  end
  if ~isempty(status.octopus)
    disp([ '  Octopus         (http://octopus-code.org/) as "' status.octopus '"' ]);
  end
  
end
% disp('Calculator executables can be specified as ''options.command=exe'' when building a model.');

%  lj (lenard-jones)
%  morse
%  eam


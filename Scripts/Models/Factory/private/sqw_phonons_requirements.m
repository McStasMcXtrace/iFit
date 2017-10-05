function status = sqw_phonons_requirements
% sqw_phonons_requirements: check for availability of ASE and MD codes
%
% MPI, EMT, GPAW, Abinit, Elk, QE, VASP
%
% returns a structure with a field for each MD software being 1 when available.
status = [];

% required to avoid Matlab to use its own libraries
if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; '; 
else           precmd=''; end

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
  [st, result] = system([ precmd status.python ' -c "from phonopy import Phonopy"' ]);
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
  for calc={'octopus','octopus_mpi'}  % octopus_mpi is obsolete, 2nd choice.
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
  
  % test for CP2K
  status.cp2k = '';
  for calc={'cp2k_shell','cp2k_shell.popt'}
    % now test executable
    [st,result]=system([ precmd 'echo EXIT | ' calc{1} ]);
    if any(st == 0)
        status.cp2k=calc{1};
        st = 0;
        break;
    end
  end
  if ~isempty(status.cp2k)
    disp([ '  CP2K            (http://www.cp2k.org/) as "' status.cp2k '"' ]);
  end
  
  % test for SIESTA
  status.siesta = '';
  for calc={'siesta'}
    % now test executable
    [st,result]=system([ precmd 'echo | ' calc{1} ]);
    if any(st == 0:1)
        status.siesta=calc{1};
        st = 0;
        break;
    end
  end
  if ~isempty(status.siesta)
    disp([ '  SIESTA          (https://departments.icmab.es/leem/siesta/) as "' status.siesta '"' ]);
  end
  
end
% disp('Calculator executables can be specified as ''options.command=exe'' when building a model.');

%  lj (lenard-jones)
%  morse
%  eam


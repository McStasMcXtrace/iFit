function status = sqw_phonons_requirements
% sqw_phonons_requirements: check for availability of ASE and MD codes
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
  status.emt='ase-run';
  disp('  EMT           only for Al,Cu,Ag,Au,Ni,Pd,Pt,H,C,N,O');
  
  % test for mpirun
  status.mpirun = '';
  for calc={'mpirun','mpiexec'}
    % now test executable
    [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
    if st == 0
        status.mpirun=calc{1};
        st = 0;
        break;
    end
  end
  
  % test for GPAW
  [st, result] = system([ precmd 'python -c "from gpaw import GPAW"' ]);
  if st == 0
    status.gpaw='gpaw-python';
    disp([ '  GPAW (http://wiki.fysik.dtu.dk/gpaw) as "' status.gpaw '"' ]);
  else
    status.gpaw='';
  end
  
  % test for NWChem
  [st, result] = system([ precmd 'python -c "from ase.calculators.nwchem import NWChem"' ]);
  status.nwchem='';
  if st == 0
    % now test executable
    % create a fake nwchem.nw
    f = tempname;
    dlmwrite([ f '.nw' ], '');
    [st,result]=system([ precmd 'nwchem ' f '.nw' ]);
    if st==0 || st==139
      status.nwchem='nwchem';
    end
    delete([ f '.*' ])
    [p,f] = fileparts(f);
    delete([ f '.db' ])
  end
  if ~isempty(status.nwchem)
    disp(['  NWChem (http://www.nwchem-sw.org/) as "' status.nwchem '"' ]);
  end
  
  % test for Jacapo
  [st, result] = system([ precmd 'python -c "from ase.calculators.jacapo import Jacapo"' ]);
  status.jacapo='';
  if st == 0
    % now test executable
    for calc={'dacapo_mpi.run','dacapo_serial.run','dacapo.run','dacapo'}
      [st,result]=system([ precmd calc{1} ]);
      if st == 0 || st == 2
          status.jacapo=calc{1};
          st = 0;
          break;
      end
    end
  end
  if ~isempty(status.jacapo)
    disp([ '  Dacapo (http://wiki.fysik.dtu.dk/dacapo) as "' status.jacapo '"' ]);
  end
  
  % test for Elk
  [st, result] = system([ precmd 'python -c "from ase.calculators.elk import ELK"' ]);
  status.elk='';
  if st == 0
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
    disp([ '  Elk (http://elk.sourceforge.net) as "' status.elk '"' ]);
  end
  
  % test for ABINIT
  [st, result] = system([ precmd 'python -c "from ase.calculators.abinit import Abinit"' ]);
  status.abinit='';
  if st == 0
    for calc={'abinit','abinis','abinip'}
      % now test executable
      [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
      if st == 0 || st == 2
          status.abinit=calc{1};
          st = 0;
          break;
      end
    end
  end
  if ~isempty(status.abinit)
    disp([ '  ABINIT (http://www.abinit.org/) as "' status.abinit '"' ]);
  end
  
  % test for QuantumEspresso
  status.quantumespresso = '';
  for calc={'pw.x','pw.exe','pw','pwscf'}
    % now test executable
    [st,result]=system([ precmd 'echo "0" | ' calc{1} ]);
    if st == 0 || st == 2
        status.quantumespresso=calc{1};
        st = 0;
        break;
    end
  end
  
  
  % test for PHON
  [st, result] = system([ precmd 'phon' ]);
  try
    delete('CRASH');
    delete('input_tmp.in');
  end
  if st == 0 || st == 2
    status.phon = 'phon';
  else
    status.phon = '';
    status.quantumespresso = '';
  end
  if ~isempty(status.quantumespresso)
    disp([ '  QuantumEspresso (http://www.quantum-espresso.org/) as "' status.quantumespresso '"' ]);
  end
  
  % test for VASP
  [st, result] = system([ precmd 'vasp' ]);
  if st == 0 || st == 2
    status.vasp = 'vasp';
  else
    status.vasp = '';
  end
  if ~isempty(status.vasp)
    disp([ '  VASP (http://www.vasp.at/) as "' status.vasp '"' ]);
  end
  % must set:
  % VASP_COMMAND=vasp
  % VASP_PP_PATH=/usr/share/vasp/pseudo
  
end
disp('Calculator executables can be specified as ''options.command=exe'' when building a model.');

%  lj (lenard-jones)
%  morse
%  eam


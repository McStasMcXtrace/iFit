function ace = read_ace(filename)
% data = read_ace(filename)
%
%   import MCNP/ACE files using PyNE
%
% (c) E.Farhi, ILL. License: EUPL.

%
% for MF7 MT4 from ACE
% lib = pyne.ace.Library('/home/farhi/Programs/MCNP/CAB-Sab/hwtr-293/hwtr-293.ace' )
% lib.read()
% scipy.io.savemat('ace.mat',lib.tables)
% matlab> ace=load('ace.mat')
%

persistent status

% ============================ Use PyNE ? ======================================
if ~exist('status') || isempty(status) || ~isstruct(status)
  status = read_ace_pyne_present;  % private function inline below
end

if status
  ace = read_ace_pyne(filename);  % private function inline below
  if ~isempty(ace), return; end
else

end

% ------------------------------------------------------------------------------
% PyNE reader routines
function status = read_ace_pyne_present
  % read_ace_pyne_present: check for availability of PyNE
  %
  % returns a flag being 1 when available.
  status = 0;

  % test for PyNE in Python
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else precmd=''; end
  [status, result] = system([ precmd 'python -c "import pyne"' ]);
  if status ~= 0  % not present
    status = 0;
    disp([ mfilename ': warning: would make good use of PyNE.' ])
  disp('  Get it at <http://pyne.io/>.');
  disp('  Package exists for Debian/Mint/Ubuntu at <http://packages.mccode.org>.');
  else
    status = 1;
  end
  
function ace = read_ace_pyne(filename)
  % read ACE file using PyNE
  
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else precmd=''; end
  
% scipy.io.savemat('ace.mat',lib.tables)
  
  % list of python 'dict' to store in the MAT file
  dict = 'atomic_relaxation, decay, fission, info, target, projectile, resonances, thermal_elastic, thermal_inelastic, reactions';
  dict = strtrim(regexp(dict, ',', 'split')); % split as words
  % create the python script to evaluate, and a lambda function to clean keys
  tmp = tempname;
  s = { 'from pyne.ace import Library ', ...
    'import scipy.io as sio ', ...
    'import re', ...
   [ 'ace=Library("' filename '") ' ], ...
    'ace.read() ', ...
    'clean = lambda varStr: re.sub("^(_)","T",  re.sub("\W|^(?=\d)","_", str(varStr))  )' };
  % each dict key is cleaned into variable name for matlab
  d = 'tables';
  s{end+1} = sprintf([ ...
    'for key in ace.%s.keys():\n' ...
    '    c=clean(key); ace.%s[c] = ace.%s.pop(key);\n' ...
    '    if type(ace.%s[c]) == dict:\n' ...
    '        for k in ace.%s[c]: ace.%s[c][clean(k)] = ace.%s[c].pop(k)' ], ...
    d,d,d,d,d,d,d);
   
  s{end+1} = sprintf('try: sio.savemat("%s_%s", ace.tables)\nexcept: None', ...
    tmp, d);

  % write script
  fid = fopen([ tmp '.py' ],'w');
  fprintf(fid, '%s\n', s{:});
  fclose(fid);
  % now evaluate in python
  try
    [status, result] = system([ precmd 'python ' tmp '.py' ]);
  catch
    status = 127; result = 'ERROR calling Python';
  end
  disp(result);
  if status ~= 0  % error occured
    disp([ mfilename ': ERROR: could not read file ' filename ' using PyNE. Trying pure Matlab method.' ]);
    ace = [];
    return
  end
  
  % import back data from MAT file
  ace = load([ tmp '_' d ]);
  ace.filename = filename;
  ace.pyne_script   = char(s);
  
  % delete temporary files created
  delete([ tmp '*' ]);
  

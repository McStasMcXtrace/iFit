function [file_out, result] = cif2hkl(varargin)
% file_out = cif2hkl(file_in, file_out, lambda, mode, verbose)
%
% Read a CIF/CFL/SHX/PCR crystallographic description
%       and generates a HKL F^2 reflection list.
%
% When the argument is a chemical formula (elements separated with spaces), or 
% a COD ID, a search in the Crystallography Open Database is made.
%
%  file_out = cif2hkl('cod: Mg O');
%  file_out = cif2hkl('cod: O3 Sr Ti');
%  file_out = cif2hkl('cod: 1532356');
%
% The formula should be given in Hill notation, e.g. C, then H, then other
% elements in alphabetical order.
%
%  file_out = cif2hkl('gui');        show a dialogue box to enter COD ID, formula, file path
%
% This requires proxy settings to be set (when behind a firewall)
%   ProxyHost='proxy.ill.fr'; % Proxy address if you are behind a proxy [e.g. myproxy.mycompany.com or empty]
%   ProxyPort=8888;           % Proxy port if you are behind a proxy [8888 or 0 or empty]
%   java.lang.System.setProperty('http.proxyHost', ProxyHost); 
%   com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(true);
%   com.mathworks.mlwidgets.html.HTMLPrefs.setProxyHost(ProxyHost);
%   java.lang.System.setProperty('http.proxyPort', num2str(ProxyPort));
%   com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPort(num2str(ProxyPort));
%
% Syntax: [file_out, result] = cif2hkl(file_in, file_out, lambda, mode, verbose)
%
% Input:
%   file_in:   name of CIF/CFL/PCR/ShellX file (char) or COD entry or COD formula
%   file_out:  name of output file, or empty for '<file_in>.laz' (char)
%   lambda:    minimum wavelength used for generation of HKL peaks. Default is 0.5  (scalar, Angs)
%   mode:      file conversion mode 'p'=powder, 'x'=single crystal, '-'=no/unactivated
%   verbose:   verbosity level, 0 or 1 (scalar)
% Output:
%   file_out:  the path to the generated output
%   result:    command result (char)
%
% Revision: $date$
% (c) E.Farhi, ILL. License: EUPL.

% we choose NOT to use the cif2hkl mex file due to allocation errors under Matlab.

persistent compiled

result = []; file_out=[];
% required to avoid Matlab to use its own libraries
exe_path = fileparts(which(mfilename));
if ismac,      precmd = [ 'DYLD_LIBRARY_PATH=' fullpath(exe_path, 'private') '; DISPLAY= ; ' ];
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; '; 
else           precmd = ''; end

% check if we use the cif2hkl executable, or need to compile
if isempty(compiled) || (nargin >0 && (strcmp(varargin{1}, 'compile') || strcmp(varargin{1}, 'check')))
  if nargin>0 && strcmp(varargin{1}, 'compile')
    compiled=cif2hkl_check_compile('compile');
  else
    compiled=cif2hkl_check_compile;
  end
end

% assemble command line
if ~isempty(compiled)
  cmd = compiled;
else
  cmd = mfilename;
end

if isempty(varargin) || strcmp(varargin{1}, 'compile') || strcmp(varargin{1}, 'check') || strcmp(varargin{1}, 'identify')
  result   = cmd;
  file_out = cmd;
  return;
end

% handle input arguments
% varargin = file_in, file_out, lambda, mode, verbose

if nargin < 1
  file_in = [];
else
  file_in = varargin{1};
end
if isempty(file_in) || ~ischar(file_in)
  error('At least one string input required. Syntax: cif2hkl(file_in, file_out, lambda, mode, verbose)')
end
if nargin < 2
  file_out = [];
else
  file_out = varargin{2};
end

if nargin < 3
  lambda = [];
else
  lambda = varargin{3};
end
if isempty(lambda) || ~isnumeric(lambda)
  lambda=0.5;
end
if nargin < 4
  mode = [];
else
  mode = varargin{4};
end
if isempty(mode) || ~ischar(mode)
  mode='p'; % powder or xtal
end
if nargin < 5
  verbose = [];
else
  verbose = varargin{5};
end
if isempty(verbose) || ~isnumeric(verbose)
  if nargout > 1, verbose=1; else verbose=0; end
end

% handle 'gui' and 'cod' entry modes
if strcmpi(file_in, 'gui') || strncmp(file_in, 'cod:', 4)
  file_in = cif2hkl_cod(file_in);
  if isempty(file_in), return; end
end

if isempty(file_out)
  file_out = [ file_in '.laz' ];
end

if ~any(strcmp(file_in, {'-h','--help','--version'}))
  
  if verbose
    cmd = [ cmd ' --verbose' ];
  end
   % cmd = [ cmd ' --lambda ' num2str(lambda) ];
  if mode(1) == '-'
    cmd = [ cmd ' --no-outout-files' ];
  else
    if mode(1) == 'p'
      cmd = [ cmd ' --powder' ];
    elseif mode(1) == 'x'
      cmd = [ cmd ' --xtal' ];
    end
    cmd = [ cmd ' --out ' file_out ];
  end
  
end

cmd = [ cmd ' ' file_in ];
disp(cmd)

% launch the command
[status, result] = system([ precmd cmd ]);
if status ~= 0
  disp(result)
  error([ mfilename ' executable not available. Compile it with: "gfortran -O2 -o cif2hkl cif2hkl.F90 -lm" in ' fullfile(fileparts(which(mfilename))) ]);
elseif any(strcmp(file_in, {'-h','--help','--version','check','compile','identify'}))
  file_out = result;
end

% ------------------------------------------------------------------------------
function compiled=cif2hkl_check_compile(compile)

  compiled = '';
  this_path = fileparts(which(mfilename));
  % required to avoid Matlab to use its own libraries
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; DISPLAY= ; '; 
  else           precmd = ''; end
  
  
  % binary external
  if ispc, ext='.exe'; else ext=''; end
  
  % look for executable, global and local
  % try in order: global(system), local, local_arch
  for try_target={fullfile(this_path, [ 'cif2hkl_' computer('arch') ext ]), ...
          fullfile(this_path, [ 'cif2hkl' ext ]), ...
          [ 'cif2hkl' ext ]}
      
          [status, result] = system([ precmd try_target{1} ]);
          if status == 0 && nargin == 0 && (~ispc || isempty(strfind(result, [ '''' try_target{1} '''' ])))
              % the executable is already there. No need to make it .
              target = try_target{1};
              compiled = target; 
              return
          end

  end
  % when we get there, target is cif2hkl_arch, not existing yet
  target = fullfile(this_path, [ 'cif2hkl_' computer('arch') ext ]);
  
  % search for a fortran compiler
  fc = '';
  for try_fc={getenv('FC'),'gfortran','g95','pgfc','ifort'}
    if ~isempty(try_fc{1})
      [status, result] = system([ precmd try_fc{1} ]);
      if status == 4 || ~isempty(strfind(result,'no input file'))
        fc = try_fc{1};
        break;
      end
    end
  end
  if isempty(fc)
    if ~ispc
      disp([ mfilename ': ERROR: FORTRAN compiler is not available from PATH:' ])
      disp(getenv('PATH'))
      disp([ mfilename ': Try again after extending the PATH with e.g.' ])
      disp('setenv(''PATH'', [getenv(''PATH'') '':/usr/local/bin'' '':/usr/bin'' '':/usr/share/bin'' ]);');
    end
    error('%s: Can''t find a valid Fortran compiler. Install any of: gfortran, g95, pgfc, ifort\n', ...
    mfilename);
  end
   
  % attempt to compile as local binary
  if isempty(dir(fullfile(this_path,mfilename))) % no executable available
    fprintf(1, '%s: compiling binary...\n', mfilename);
    cmd = {fc, '-o', target, ...
       fullfile(this_path,'cif2hkl.F90'), '-lm', '-ffree-line-length-0'}; 
    disp([ sprintf('%s ', cmd{:}) ]);
    [status, result] = system([ precmd sprintf('%s ', cmd{:}) ]);
    if status ~= 0 % not OK, compilation failed
      disp(result)
      warning('%s: Can''t compile cif2hkl.F90 as binary\n       in %s\n', ...
        mfilename, fullfile(this_path));
    else
      delete(fullfile(this_path, '*.mod'));
      delete('*.mod');
      compiled = target;
    end
  end
  
% gui and cod entries ----------------------------------------------------------
function file = cif2hkl_cod(file)
  
  if strcmpi(file, 'gui')
    % pops-up a dialogue to enter a chemical formula
    prompt={'Enter a {\color{red}chemical formula} in Hill notation (C, then H, then other elements in alpha order, separated with spaces) or {\color{red}COD ID} to search on {\color{blue}crystallography.net} or a {\color{red}file name}:'};
    name=[ mfilename ': Input chemical formula/COD ID' ];
    numlines=1;
    defaultanswer={'O3 Sr Ti'};
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter = 'tex';

    file=inputdlg(prompt,name,numlines,defaultanswer,options);
    if isempty(file), return; else file = [ 'cod: ' char(file) ]; end
  end
  
  
  if ~strncmp(file, 'cod:', 4)
    return
  else
    file = file(5:end);
  end
  
  file = strtrim(file);
  
  % is this a COD ID ?
  if ~isnan(str2double(file))
    cod = [ file '.cif' ];
  else
    formula = strrep(strtrim(file), ' ', '%20');  % change spaces from formula to cope with COD query
    
    % query COD
    disp([ mfilename ': querying COD at http://www.crystallography.net/cod/result.php?formula=' formula ]);
    try
      cod   = urlread([ 'http://www.crystallography.net/cod/result.php?formula=' formula ]);
    catch
      disp([ mfilename ': It seems I can not reach www.crystallography.net !' ]);
      disp('>>> if you are behind a Proxy, you MUST set from the Matlab/iFit prompt e.g.:');
      disp('  ProxyHost=''proxy.ill.fr'' % no need for "http://" there');
      disp('  ProxyPort=8888');
      disp('  java.lang.System.setProperty(''http.proxyHost'', ProxyHost); ')
      disp('  com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(true);');
      disp('  com.mathworks.mlwidgets.html.HTMLPrefs.setProxyHost(ProxyHost);');
      disp('  java.lang.System.setProperty(''http.proxyPort'', num2str(ProxyPort));');
      disp('  com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPort(num2str(ProxyPort));');
      error([ mfilename ': Network seems unreachable. Check connection or Proxy settings.' ]);
    end
    cod   = textscan(cod, '%s','Delimiter',sprintf('\n\r'));
    cod   = cod{1};
    cod   = cod(find(~cellfun(@isempty, strfind(cod, 'CIF'))));
    % have to read lines until '</tr>'
    index = strfind(cod, '</tr>');
    cod   = cod(find(~cellfun(@isempty, index)));
    index = cell2mat(index);
    for l=1:numel(cod)
      this = cod{l};
      this = this(1:(index(l)-1));
      % remove some of the links '<a href="result.php
      i1 = strfind(this, '<a href="result.php?spacegroup');
      i2 = strfind(this, '>');
      if ~isempty(i1) && ~isempty(i2)
        i1 = i1(1); i2=i2(find(i2 > i1, 1));
        this(i1:(i2)) = [];
      end
      i3 = strfind(this, '<a href="result.php?journal');
      if ~isempty(i3)
        i3 = i3(1);
        this(i3:end) = [];
      end
      this = strrep(this, '<a href="', '');
      % clean up each CIF line: remove <td> </td> <br/> <i> </i> <b> </b> ...
      for tok={'<td>', '</td>', '<br/>', '<i>', '</i>', '<b>', '</b>', '%20','</a>', '<a href="', '">'}
        this = strrep(this, tok{1}, ' ');
      end
      % the first token in each line is now the COD number
      cod{l} = this;
    end
    % pop-up  dialogue to choose when more than one entry
    if isempty(cod), file=''; return; end
    if numel(cod) > 1
      selection = listdlg('ListString', cod, 'ListSize', [ 300 160 ], ...
        'Name', [ mfilename ': Crystallography Open Database entries for ' file ], ...
        'SelectionMode', 'single', ...
        'PromptString', { [ 'Here are the entries for ' file ]; ...
        'from the Crystallography Open Database (COD) <http://www.crystallography.net>.'; ...
        'Each entry shows the COD ID, spacegroup, cell parameters (a,b,c,alpha,beta,gamma) and title.'; ...
        'Please choose one of them.' });
      if isempty(selection),    file=''; return; end % cancel
    else selection = 1; end
    if numel(cod) && iscellstr(cod)
      cod = cod{selection};
    end
  end
  
  % get the COD entry (first token from 'cod' string)
  cod_id = strtok(cod);
  disp([ mfilename ': getting http://www.crystallography.net/cod/' cod_id ]);
  file = urlread([ 'http://www.crystallography.net/cod/' cod_id ]);
  % copy that file locally in temp dir
  try
    d = tempname;
    mkdir(d);
    fid=fopen(fullfile(d, cod_id), 'w');
    if fid==-1, disp([ mfilename ': could not write ' cod_id ]); end
    fwrite(fid, file); fclose(fid);
    file = fullfile(d, cod_id);
  end

  

function fhandle = ifitmakefunc(fun, descr, pars, expr, defPars, constraint)
% fhandle = ifitmakefunc(fun, descr, pars, expr, defpars, constraint) : build a fit function/model
%
%   iFit/ifitmakefunc fit function/model builder.
%     when input parameters are missing, a dialog pops-up to request the information
%       all arguments are optional (can be left empty), except the expression ;
%     when only the first argument is specified as an expression (using p and x,y,z,t)
%       a model is built from the expression analysis.
%         ifitmakefunc(EXPR)
%     when only the first two arguments are specified as an expression (using p 
%       and x,y,z,t) and a constraints a model is built from the expression analysis.
%         ifitmakefunc(EXPR, CONSTRAINT)
%     when only the first argument is specified as a structure, it should define fields
%       fun.function:    the function name (used for definition and file storage)
%       fun.Description: the descriptionnof the function
%       fun.Parameters:  the parameter names as words separated with spaces
%       fun.Guess:       the default parameters, or 'automatic'
%       fun.Expression:  the expression of the function value
%       fun.Constraint:  the expression executed before the function evaluation
%     to make the new function permanent, copy it into [ ifitpath '/iFuncs' ]
%     The list of all available function can be obtained with the 'fits(iData)' 
%
% input:  FUN:     function name or expression (single word or expression or structure)
%         DESCR:   description of the function (string)
%         PARS:    name of parameters, as a single string of words 'a b c ...'
%         EXPR:    expression of the function value, using 'p' vector as parameter values
%         DEFPARS: default parameter values (guess, optional, leave empty for automatic guess)
%         CONSTRAINT: expression to execute before the function evaluation
%           which may contain any parameter constraints such as 'p(1)=p(5)'
%
% output: fhandle: function handle to the new function, which is also stored locally
% 
% Version: $Revision: 1.7 $
% See also iData, gauss

fhandle = [];
NL      = sprintf('\n');

if nargin == 1 && strcmp(fun,'identify')
  y.Type           = 'iFit function builder';
  y.Name           = [ 'iFit function maker [' mfilename ']' ];
  y.Parameters     = '';
  y.Dimension      = 0;         % dimensionality of input space (axes) and result
  y.Guess          = [];        % default parameters
  y.Axes           = {};        % the axes used to get the values
  y.Values         = [];        % default model values=f(p)
  y.function       = mfilename;
  fhandle=y;
  return
elseif nargin == 1 && strcmp(fun,'plot')
  text(0.1, 0.5, 'Function maker');
  title(mfilename);
  return
elseif nargin == 1 && ~isvarname(fun)
  % special case when only the expression is given
  expr   = fun;
  fun    = '';
  fun    = '';
  descr  = '';
  pars   = '';
  defPars= '';
  constraint='';
elseif nargin == 2 && ~isvarname(fun)
  % special case when only the expression and the constraint are given
  expr   = fun;
  constraint=descr;
  fun    = '';
  fun    = '';
  descr  = '';
  pars   = '';
  defPars= '';
else
  % general case: handle incomplete input and pops-up dialog
  if nargin == 1 && isstruct(fun)
    if isfield(fun, 'defPars'),         defPars=fun.defPars;
    elseif isfield(fun, 'x0'),          defPars=fun.x0;
    elseif isfield(fun, 'Guess'),       defPars=fun.Guess; end
    if isfield(fun, 'Expression'),      expr=fun.Expression; end
    if isfield(fun, 'pars'),            pars=fun.pars;
    elseif isfield(fun, 'Parameters'),  pars=fun.Parameters; end
    if isfield(fun, 'Name'),            descr=fun.Name; 
    elseif isfield(fun, 'Description'), descr=fun.Description; end
    if isfield(fun, 'Constraint'),      constraint=fun.Constraint; end
    if isfield(fun, 'function'),        fun=fun.function; fun=[]; fun=tmp; end
  else
    if nargin < 6,  constraint=''; end
    if nargin < 5,  defPars='automatic'; end
    if nargin < 4,  expr = 'p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4)'; end
    if nargin < 3   pars = 'Amplitude Centre HalfWidth Background'; end
    if nargin < 2,  descr= 'Gaussian function'; end
    if nargin < 1   fun  = 'gauss'; end
  end
  
  if nargin < 4
    % request input dialog
    prompt    = { [ '{\bf Function name}' NL '(a single word, optional)' ], ...
                  [ '{\bf Description of the fit model}' NL '(a character string, optional)' ], ...
                  [ '{\bf Model parameters names}' NL '(single names separated by spaces, optional)' ], ...
                  [ '{\bf Value of the function {\color{red} required}}' NL '(expression using parameters from vector {\color{blue} p(1), p(2)}, ... and axes {\color{blue} x, y, z, t}, ...)' ], ...
                  [ '{\bf Default parameter values}' NL '(vector,  e.g [1 2 ...], leave empty to use automatic guess)' ], ...
                  [ '{\bf Constraint}' NL '(any expresion executed before the function Value, optional)' ]};
    dlg_title = 'iFit: Make fit function';
    num_lines = [ 1 1 1 3 1 3]';
    defAns    = {fun, descr, pars, expr, defPars, constraint};
    options.Resize      = 'on';
    options.WindowStyle = 'normal';   
    options.Interpreter = 'tex';
    answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);
    if isempty(answer), 
      return; 
    end
    % extract results
    fun   = answer{1};
    descr = answer{2};
    pars  = answer{3};
    expr  = answer{4};
    defPars=answer{5};
    constraint=answer{6};
  end % nargin < 4
end % else 

% default function name, and possibly extract path for storage
if isempty(fun)
  fun = clock;
  fun = sprintf('fun_%i', fix(fun(6)*1e4));
end
[ fun_path, fun ]   = fileparts(fun);
fun                 = genvarname(fun);
nb_pars             = 0;
dim                 = 0;

% default parameter names
if findstr(defPars, 'auto'), defPars = ''; end

if ~isempty(expr)
  % compute the number of parameters and corresponding parameter names
  nb_pars = findstr(expr, 'p('); % guess may be wrong when matching e.g. 'exp('...
  % we look if this is really a p(n) syntax
  n = [];
  for index=1:length(nb_pars)
    if nb_pars(index)>2 && isstrprop(expr(nb_pars(index)-1),'alpha'), 
      continue; 
    end
    token          = expr((nb_pars(index)+2):end);
    closing_parent = strfind(token, ')');
    if isempty(closing_parent), closing_parent=length(token)+1; end
    try
      s = eval(token(1:(closing_parent-1)));
    catch
      s=[];
    end
    if length(s) >= 1 ... 
      n=[n s]; 
    end
  end
  if ~isempty(n), nb_pars = max(n);
  else nb_pars=0; end
  if ~isempty(pars)
    % check if the number of parameters used in the expression matches the parameter names
    if length(strread(pars, '%s','delimiter',' ;,''{}')) ~= nb_pars
      fprintf(1, [ '%s: Warning: iFit model ''%s''' NL ...
                   '    ''%s''' NL ...
                   '  The number of parameters %d used in the expression' NL ...
                   '  may not be the same as the number of parameter names %d\n' ], ...
                   mfilename, fun, expr, nb_pars, length(strread(pars, '%s','delimiter',' ;,''{}')));
      if length(strread(pars, '%s','delimiter',' ;,''{}')) < nb_pars
        pars = []; % will guess parameter names
      end
    end
  end
  
  if ~isempty(defPars)
    % check if the number of parameters used in the expression matches the parameter default values
    if length(defPars) ~= nb_pars
      fprintf(1, [ '%s: Warning: iFit model ''%s''' NL ...
                   '    ''%s''' NL ...
                   '  The number of parameters %d used in the expression' NL ...
                   '  may not be the same as the number of default parameter values %d\n' ], ...
                   mfilename, fun, expr, nb_pars, length(defPars));
    end
  end
  
  % compute model dimensionality from occurence with x,y,z,t
  if length(findstr(expr, 'x'))
    dim = 1;
    if length(findstr(expr, 'y'))
      dim = 2;
      if length(findstr(expr, 'z'))
        dim = 3;
        if length(findstr(expr, 't'))
          dim = 4;
        end
      end
    end
  end % if 'x' (1d)
end % expr

% return when model can not be defined
if isempty(expr) || dim == 0 || nb_pars == 0
  return
end
% dim:     holds model dimensionality
% nb_pars: holds number of parameter used in expression

% default description
if isempty(descr), descr=expr; end

% default parameter names
% try to be clever by an analysis of the expression...
if isempty(pars)
  for index=1:nb_pars
    namp = [ 'p(' num2str(index) ')' ];
    if ~isempty(findstr(expr, [ '*' namp ])) || ~isempty(findstr(expr, [ namp '*' ]))
      name = [ 'Amplitude_' num2str(index) ];
    elseif ~isempty(findstr(expr, [ '-' namp ]))
      name = [ 'Centre_' num2str(index) ];
    elseif ~isempty(findstr(expr, [ '+' namp ])) || ~isempty(findstr(expr, [ namp '+' ]))
      name = [ 'Constant_' num2str(index) ];
    elseif ~isempty(findstr(expr, [ '/' namp ])) || ~isempty(findstr(expr, [ namp '^2' ]))  
      name = [ 'Width_' num2str(index) ];
    else
      name = [ fun '_p' num2str(index) ];
    end
    if index==1, 
      pars = [ '''' name '''' ];
    else
      pars = [ pars ', ''' name '''' ];
    end
  end
else
  pars     = strread(pars, '%s','delimiter',' ;,''{}'); % is a cellstr
  % re-arrange the parameter names as a string prepresenting a cellstr
  new_pars = '';
  for index=1:length(pars)
    if ~isempty(pars{index})
      new_pars = [ new_pars ' ''' genvarname(pars{index}) '''' ];
    end
  end
  pars = new_pars;
end

% ============================ CREATE model ====================================

% load template from [ fileparts(which(mfilename)) filesep 'private' filespe 'template.txt' ]
template = fullfile(fileparts(which(mfilename)), 'private', 'template.txt');
[fid, message] = fopen(template);
if fid == -1
  error(sprintf('Can not open %s for reading\n  %s\n', template, message));
end

% read file contents as a string
template = fread(fid, Inf, 'uint8=>char')';
fclose(fid);

% replace symbols:
% $FUN:   function name (single word for the file name)
% $DESCR: description of the function
% $PARS:  name of parameters, as 'a','b', ...
% $EXPR:  expression of the function value, using 'p' vector as parameter values
template = strrep(template, '$FUN',   fun);
template = strrep(template, '$DESCR', descr);
template = strrep(template, '$PARS',  pars);
if ~isempty(constraint)
  template = strrep(template, '$EXPR',  [ constraint ';' NL '  signal = ' expr ]);
else
  template = strrep(template, '$EXPR',  [ '  signal = ' expr ]);
end
whos dim
template = strrep(template, '$DIM',   num2str(dim));
ax = '';
if dim>=1, ax = [ ax   'x' ]; end
if dim>=2, ax = [ ax ', y' ]; end
if dim>=3, ax = [ ax ', z' ]; end
if dim>=4, ax = [ ax ', t' ]; end
template = strrep(template, '$AXES',  ax);
% handle default parameters: static or automatic
if ~isempty(defPars)
  template = strrep(template, '$DEFPARS', num2str(defPars));
else
  template = strrep(template, '$DEFPARS', 'iFuncs_private_guess(x(:), signal(:), info.Parameters)');
end


% add the private functions
% iFit/iFuncs/private/iFuncs_private_findpeaks.m
% iFit/iFuncs/private/iFuncs_private_guess.m
line = '% ==============================================================================';
tmp = fullfile(fileparts(which(mfilename)), 'private', 'iFuncs_private_findpeaks.m');
fid = fopen(tmp);
tmp = fread(fid, Inf, 'uint8=>char')';
fclose(fid);
template = [ template NL tmp ];

tmp = fullfile(fileparts(which(mfilename)), 'private', 'iFuncs_private_guess.m');
fid = fopen(tmp);
tmp = fread(fid, Inf, 'uint8=>char')';
fclose(fid);
template = [ template NL line NL tmp ];


% create the function: write it
[fid, message] = fopen([ fun '.m' ], 'w');
if fid == -1
  error(sprintf('Can not create function %s for writing\n  %s\n', fullfile(fun_path, [ fun '.m' ]), message));
end

fwrite(fid, template, 'char');
fclose(fid);

% display information and force to rehash/register the function
fprintf(1, [ '%s: Wrote function signal=%s(p, %s)\n'...
  '%% %s\n%% %d parameter(s): %s\n' ], ...
  mfilename, fun, ax, descr, nb_pars, ...
  pars);
if ~isempty(constraint)
  fprintf(1, '  %s;\n', constraint);
end
% display information and force to rehash/register the function
fprintf(1, '  signal=%s;\n\nStored in: %s\n', expr, which([ fun '.m' ]));

% create the handle
fhandle = str2func(fun);


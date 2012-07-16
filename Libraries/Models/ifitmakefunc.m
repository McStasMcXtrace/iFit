function f = ifitmakefunc(fun, descr, pars, expr, defPars, constraint)
% f = ifitmakefunc(fun, descr, pars, expr, defpars, constraint) : build a fit function/model
%
%   iFit/ifitmakefunc fit function/model builder.
%     when input parameters are missing, a dialog pops-up to request the information
%       all arguments are optional (can be left empty), except the expression ;
%     when only the first argument is specified as an expression (using p and x,y,z,t)
%       a model is built from the expression analysis.
%         ifitmakefunc(EXPR)
%     when only the first two arguments are specified as an expression (using p 
%       and x,y,z,t) and a constraints, a model is built from the expression analysis.
%         ifitmakefunc(EXPR, CONSTRAINT)
%     when only the first argument is specified as a structure, it should define fields
%       fun.function:    the function name (used for definition and file storage)
%       fun.Description: the description of the function
%       fun.Parameters:  the parameter names as words separated with spaces
%       fun.Guess:       the default parameters, or 'automatic'
%       fun.Expression:  the expression of the function value
%       fun.Constraint:  the expression executed before the function evaluation
%     to make the new function permanent, copy it into [ ifitpath '/Models' ]
%     The list of all available function can be obtained with the 'fits(iData)' 
%
% WARNING: the 1D and 2D cases work fine, but higher dimensions (using z,t...) are NOT validated.
%
% input:  FUN:     function name or expression (single word or expression or structure)
%         DESCR:   description of the function (string)
%         PARS:    name of parameters, as a single string of words 'a b c ...'
%         EXPR:    expression of the function value, using 'p' vector as parameter values
%         DEFPARS: default parameter values (guess, optional, leave empty for automatic guess)
%         CONSTRAINT: expression to execute before the function evaluation
%           which may contain any parameter constraints such as 'p(1)=p(5)'
%
% output: f: new model object
% 
% Version: $Revision: 1.12 $
% See also iData, gauss, iFunc

f       = [];
NL      = sprintf('\n');

if nargin == 1 && strcmp(fun, 'identify')
  return
elseif nargin == 1 && ~isvarname(fun) && (ischar(fun) || isa(fun,'function_handle'))
  % special case when only the expression is given
  expr   = fun;
  fun    = '';
  descr  = '';
  pars   = '';
  defPars= '';
  constraint='';
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
    if nargin < 4,  expr = '@(p,x)p(1)*exp(- 0.5*((x-p(2))/p(3)).^2) + p(4)'; end
    if nargin < 3   pars = 'Amplitude Centre HalfWidth Background'; end
    if nargin < 2,  descr= 'Gaussian function'; end
    if nargin < 1   fun  = 'gauss'; end
  end
  
  if nargin < 4
    % request input dialog
    prompt    = { [ '{\bf Function name}' NL '(a single word, optional)' ], ...
                  [ '{\bf Description of the fit model}' NL '(a character string, optional)' ], ...
                  [ '{\bf Model parameters names}' NL '(single names separated by spaces, optional)' ], ...
                  [ '{\bf Value of the function {\color{red} required}}' NL '(expression using parameters from vector {\color{blue} p(1), p(2)}, ... and axes {\color{blue} x, y, z, t}, ...). It may be a function handle {\color{blue} @(p,x,y...)expression}.' ], ...
                  [ '{\bf Default parameter values}' NL '(vector, expression, function handle {\color{blue} @(x,..signal)guess},  e.g [1 2 ...], leave empty for automatic guess)' ], ...
                  [ '{\bf Constraint}' NL '(any expresion or function handle {\color{blue} @(p,x,...) constraint} executed before the function Value, optional)' ]};
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

% check expression, guess and constraint for function_handle
try
  if isa(eval(expr), 'function_handle')
    expr = eval(expr);
  end
end
try
  if isa(eval(constraint), 'function_handle')
    constraint = eval(constraint);
  end
end
try
  if isa(eval(defPars), 'function_handle')
    defPars = eval(defPars);
  end
end
% default parameter names
if strncmp(defPars, 'auto',4), defPars = ''; end
if ischar(pars) 
  pars     = strread(pars, '%s','delimiter',' ;,''{}'); % is a cellstr
end
pars=pars(:)';

% create a structure before we build the object
f.Name        = fun;
f.Guess       = defPars;
f.Expression  = expr;
f.Parameters  = pars;
f.Constraint  = constraint;
f.Description = descr;

f = iFunc(f);


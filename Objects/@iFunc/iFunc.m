function a = iFunc(varargin)
% model = iFunc(...): create a fit model function
%
% Any Fit Model function can be created from a mathematical expression, a
% structure, a function handle or an other object.
%
% The model can store the following information members:
%   Expression:       The expression for the function (string) or function handle
%                       signal=Expression(p, x,y, ...)
%   Guess:            Vector, expression or function to evaluate in order to obtain
%                       guessed parameters from axes x,y,z, ... and signal
%                       p=Guess(x,y,z, ..., signal) should return a vector. NaN values 
%                       are not set.
%   Name:             A short name for the function
%   ParameterValues:  Some parameter values to be used e.g. as starting parameters
%                       for a fit procedure or a plot
%   Description:      A string describing the model
%   Parameters:       The Parameter names (cellstr)
%   Constraint:       An expression or a function handle with syntax 
%                       p=Constraint(p, x,y,...). NaN values are not changed.
%
% From a character string
%   the Expression should make use of x,y,z,t,u to denote axes of rank 1-5,
%   and the model Parameters are specified using 'p(n)' vector elements.
%
% From a structure, with iFunc object fields (see above) and alias fields:
%   x0          -> Guess
%   objective   -> Expression
%   pars        -> Parameters
%   function    -> Name
%
% From a function handle
%  The function should have syntax model(p, x,y,...) and return the model value.
% 
% input:  s: iFunc, string, structure, function handle
% output: b: model object (iFunc)
% ex:     b=iFunc(@(p,x) p(1)*x+p(2)); 
%         b=iFunc('p(1)*x+p(2)');
%         b=iFunc('signal=p(1)*x+p(2);');
%
% Version: $Revision: 1.3 $
% See also iFunc, iFunc/feval, iFunc/plot, iFunc/fit


persistent id % the id number for all new objects

a = [];

if nargin == 0  % create the empty iFunc object structure
  % create a new iFunc object
  a.Tag         = 0;
  a.Date        = clock;  % creation date
  a.Name        = ''; % the function Name
  a.Description = ''; % the Description
  a.Parameters  = {}; % the function parameters, a cellstr or words or a single
                      % string of words separated by spaces, or a structure with parameter values
  a.Guess       = 'automatic'; % the default parameter set or expression or function handle
                      % char or function_handle(x,y,..., signal, Parameters{})
                      % or 'automatic'
  a.Expression  = ''; % the expression to evaluate to get the function value
                      % char or function_handle(p, x,y, ...)
  a.Constraint  = ''; % code to evaluate before computing the Expression
  a.Dimension   = 0;  % function dimensionality (1,2,3,4...) 0=scalar=empty
  a.ParameterValues = [];
  a.Eval        = ''; % code to evaluate for the model value
  a.UserData    = '';

  if isempty(id),   id=0; end
  if id > 1e6, id=0; end
  if id <=0, 
    id = a.Date;
    id = fix(id(6)*1e4); 
  else 
    id=id+1;
  end

  a.Tag      = [ 'iF' sprintf('%0.f', id) ];
  a = class(a, 'iFunc');
elseif length(varargin) > 1 % import data to create the object array
  % do we have a name/value set ?
  if rem(length(varargin),2) == 0 && iscellstr(varargin(1:2:(end-1)))
    try
      a = iFunc(struct(varargin{:}));
      return
    end
  end
  
  % can be a structure, or an expression, or a function_handle
  for index=1:length(varargin)
    a = [ a iFunc(varargin{index}) ];
  end
  return
else   % import data to create a single object
  % can be a structure, or an expression, or a function_handle
  this = varargin{1};
  
  if ischar(this) % ----------------------------------------------------------
    % check if this is a predefined function
    if exist(this) == 2 
      a=iFunc(feval(this));
      a.Name = this;
    else
      % from a string/expression: analyse the expression to get the parameter names,
      % dimension and parameter names
      a=iFunc;
      a.Expression = this;
    end
  
  elseif isa(this, 'iFunc') % ------------------------------------------------
    % make a copy of the initial object
    a = this;
  elseif isstruct(this) % ----------------------------------------------------
    a = iFunc;
    % identify if some of the structure fields can be directly used, and overlayed
    if isfield(this, 'x0'),              a.Guess=this.x0;
    elseif isfield(this, 'Guess'),       a.Guess=this.Guess; end
    if isfield(this, 'objective'),       a.Expression=this.objective; end
    if isfield(this, 'Expression'),      a.Expression=this.Expression; end
    if isfield(this, 'pars'),            a.Parameters=this.pars;
    elseif isfield(this, 'Parameters'),  a.Parameters=this.Parameters; end
    if isfield(this, 'function'),        a.Name=this.function; 
    elseif isfield(this, 'Name'),        a.Name=this.Name; end
    if isfield(this, 'Description'),     a.Description=this.Description; end
    if isfield(this, 'Constraint'),      a.Constraint=this.Constraint; end
    
  elseif isa(this, 'function_handle') % --------------------------------------
    a = iFunc;  % empty object
    if abs(nargout(this)) < 1
      error(['iFunc:' mfilename ], '%s: function %s should return at least one output value.\n  signal=f(p, axes{:}, additional_arguments{:})', ...
        mfilename, func2str(this));
    end
    if ~abs(nargin(this)) || ~abs(nargout(this))
      error(['iFunc:' mfilename ], '%s: function %s should use at least one input and output arguments.\n  signal=f(p, x,y,z, ...)', ...
        mfilename, func2str(this));
    end
    if exist(func2str(this))
      a = feval(this);
    else
      a.Expression  = this;
      if a.Dimension == 0
        a.Dimension   = nargin(this) - 1;
      end
    end
  else
    error(['iFunc:' mfilename ], [mfilename ': import of ' inputname(1) ' of class ' class(this) ' is not supported. Use struct, function handle, char, iFunc object.' ]);
  end
  
  % check parameter names wrt expression and dimensionality, ...
  a = iFunc_private_check(a);
  
end % if

% ------------------------------------------------------------------------------
% iFunc_private_checkexpr: check the function expression/object
function a = iFunc_private_check(a)

  nb_pars             = 0;
  dim                 = 0;
  
  % Expression can be of char, cellstr and function_handle
  if isa(a.Expression,'function_handle')
    expr = func2str(a.Expression);
  elseif ischar(a.Expression) || iscellstr(a.Expression) 
    expr = char(a.Expression);
  else
    error(['iFunc:' mfilename ], [mfilename ': the model Expression should be a char or function_handle or cellstr, not a class ' ...
      class(a.Expression) '.' ]);
  end
  % handle multiple lines in Expression -> single char with \n at end of lines (except last)
  n_expr = '';
  for index=1:size(expr, 1)
    if index == size(expr, 1)
      n_expr = [ n_expr deblank(expr(index,:)) ';' ];
    else
      n_expr = [ n_expr deblank(expr(index,:)) sprintf(';\n') ];
    end
  end
  expr = n_expr;
  % Constraint can be of char, cellstr
  if ischar(a.Constraint) || iscellstr(a.Constraint) 
    const = char(a.Constraint);
  else
    error(['iFunc:' mfilename ], [mfilename ': the model Constraint should be a char or cellstr, not a class ' ...
      class(a.Constraint) '.' ]);
  end
  % handle multiple lines in Expression -> single char with \n at end of lines (except last)
  n_const = '';
  for index=1:size(const, 1)
    if index == size(const, 1)
      n_const = [ n_const deblank(const(index,:)) ';' ];
    else
      n_const = [ n_const deblank(const(index,:)) sprintf(';\n') ];
    end
  end
  const = n_const;
  
  if ischar(a.Parameters)
    pars              = strread(a.Parameters, '%s','delimiter',' ;,''"{}'); % is a cellstr
  elseif isstruct(a.Parameters)
    pars              = fieldnames(pars);
  elseif iscellstr(a.Parameters)
    pars              = a.Parameters;
  end
  
  if ~isempty(pars) && ~iscellstr(pars)
    disp([mfilename ': the model parameters should be a char or structure or cellstr, not a class ' ...
      class(pars) '. Setting new parameter names.' ]);
    pars = {};
  end

  if strncmp(a.Guess, 'auto',4), a.Guess = ''; end

  if ~isempty(expr)
    % first count the initial p(n) fields and get the used indices
    % compute the number of parameters and corresponding parameter names
    % regexp: p(...) with '[0-9:]' characters inside
    nb_pars = regexp(expr,'\<p\(([\[0-9\:\]]+)\)','tokens'); % return the tokens
    if ~isempty(nb_pars)
      % assemble all matches
      n = '';
      for index=1:length(nb_pars)
        this = nb_pars{index};
        n = [ n this{1}  ' ' ];
      end
      nb_pars = unique(str2num([ '[' n ']' ]));
      used_pars=zeros(1,max(nb_pars));
      used_pars(nb_pars) = 1;
    else
      used_pars=[];
      nb_pars=[];
    end

    if ~isempty(nb_pars), 
      nb_pars = max(nb_pars);
    else nb_pars=0; end
    NL = sprintf('\n');
    if ~isempty(pars)
      % check if the number of parameters used in the expression matches the parameter names
      if length(pars) ~= nb_pars
        fprintf(1, [ '%s: Warning: iFit model ''%s''' NL ...
                     '    ''%s''' NL ...
                     '  The number of parameters %d used in the expression is not' NL ...
                     '  the same as the number of parameter names %d. Fixing.\n' ], ...
                     mfilename, a.Tag, expr, nb_pars, length(pars));

        pars = []; % will guess parameter names
      end
    end
    
    if ~isempty(a.Guess) && ~isa(a.Guess, 'function_handle')
      % check if the number of parameters used in the expression matches the parameter default values
      if ~ischar(a.Guess) && length(a.Guess) ~= nb_pars
        fprintf(1, [ '%s: Warning: iFit model ''%s''' NL ...
                     '    ''%s''' NL ...
                     '  The number of parameters %d used in the expression' NL ...
                     '  may not be the same as the number of default parameter values %d\n' ], ...
                     mfilename, a.Tag, expr, nb_pars, length(a.Guess));
      end
    end
    
    % compute model dimensionality from occurence with x,y,z,t
    if length(regexp(expr, '\<x\>'))
      dim = 1;
      if length(regexp(expr, '\<y\>'))
        dim = 2;
        if length(regexp(expr, '\<z\>'))
          dim = 3;
          if length(regexp(expr, '\<t\>'))
            dim = 4;
            if length(regexp(expr, '\<u\>'))
              dim = 5;
            end
          end
        end
      end
    end % if 'x' (1d)
  end % if ‪~empty(expr)

  % return when model can not be defined
  if isempty(expr) || nb_pars == 0
    warning('iFunc:emptyModel', 'Expression does not contain any axes (x,y,z...) or parameters (p). Empty model.'); 
    return
  end
  % dim:     holds model dimensionality
  % nb_pars: holds number of parameter used in expression
  if ~isa(a.Expression, 'function_handle') 
    % check that Expression can define the 'signal'
    e = textscan(expr,'%s','Delimiter',sprintf('\n\r\f'),'MultipleDelimsAsOne',1); e=e{1};
    has_signal = 0;
    for index=1:length(e)
      e{index} = strtrim(e{index});
      if ~isempty(regexp(e{index}, '\<signal\>\s*=')), has_signal = 1; end
      if index == length(e) && ~has_signal
        e{index} = sprintf('signal = %s', e{index});
        has_signal = 1;
      end
    end
    e = strrep(e, ';;','; ');
    a.Expression = char(e);
  end % else keep function handle value
  a.Dimension  = dim;

  % default parameter names
  % try to be clever by an analysis of the expression...
  amp=0; cen=0; bkg=0; wid=0;
  if isempty(pars)
    pars = {};
    for index=1:nb_pars
      namp = [ 'p(' num2str(index) ')' ];
      if ~isempty(findstr(expr, [ '*' namp ])) || ~isempty(findstr(expr, [ namp '*' ]))
        if ~amp
          name = 'Amplitude'; amp=1;
        else
          name = [ 'Amplitude_' num2str(index) ];
        end
      elseif ~isempty(findstr(expr, [ '-' namp ]))
        if ~cen
          name = 'Centre'; cen=0;
        else
          name = [ 'Centre_' num2str(index) ];
        end
      elseif ~isempty(findstr(expr, [ '+' namp ])) || ~isempty(findstr(expr, [ namp '+' ]))
        if ~bkg
          name = 'Constant'; bkg=1;
        else
          name = [ 'Constant_' num2str(index) ];
        end
      elseif ~isempty(findstr(expr, [ '/' namp ])) || ~isempty(findstr(expr, [ namp '^2' ]))  
        if ~wid
          name = 'Width'; wid=1;
        else
          name = [ 'Width_' num2str(index) ];
        end
      else
        name = [ a.Tag '_p' num2str(index) ];
      end
      pars{index} = name;
    end
  else
    % re-arrange the parameter names as a string representing a cellstr
    for index=1:length(pars)
      if ~isempty(pars{index})
        u = pars{index};
        u(~isstrprop(u,'print'))=' ';
        pars{index} = u;
      end
    end
  end
  
  % check parameter values (if any defined yet): should be as many as parameter names
  val = [];
  if isstruct(a.Parameters) && length(fieldnames(a.Parameters)) == length(pars)
    val = cell2mat(struct2cell(a.Parameters));
  elseif ~isempty(a.ParameterValues)
    val = a.ParameterValues;
  end
  a.Parameters      = pars;
  if ~isempty(val)
    if length(pars) ~= length(val)
      fprintf(1,'%s: model %s: the number of model parameter values %d does not match the number of parameter names %d.\n', mfilename, a.Tag, length(val), length(pars));
      val = zeros(size(pars));
    end
    a.ParameterValues = val;
  end
  
  if isempty(a.Eval)
    a.Eval = char(a);
  end
  
  % create a default Name for Labels
  if isempty(a.Name)
    e = textwrap(cellstr(char(a.Expression)),80);
    if length(e) > 3, e=e(1:3); end
    e = sprintf('%s', e{:}); if e(end) == ';', e(end)=''; end
    token = regexp(e, '\<signal\>\s*=','match');
    if ~isempty(token), e = strrep(e, token{1}, ''); end
    if length(e) > 20, e=[ e(1:17) '...' ]; end
    a.Name = e;
  end
  
  % handle non char fields to convert to single strings
  if ~isvector(a.Description), a.Description=char(a.Description)'; a.Description = a.Description(1:end); end
  if ~isvector(a.Expression),  a.Expression =char(a.Expression)';  a.Expression  = a.Expression(1:end); end
  if ~isvector(a.Constraint),  a.Constraint =char(a.Constraint)';  a.Constraint  = a.Constraint(1:end); end
  if ~isvector(a.Guess),       a.Guess      =char(a.Guess)';       a.Guess       = a.Guess(1:end); end

% end iFunc_private_check

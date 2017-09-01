classdef iFunc
% model = iFunc(...): create a fit model function
%
% Any Fit Model function can be created from a mathematical expression, a
% structure, a function handle or an other object.
%
% The model can store the following information members:
%   Expression:       The expression for the function (string) or function handle
%                       signal=Expression(p, x,y, ...)
%                         The parameter names surrounded by "" are replaced
%                         by the corresponing p(n)
%                       the iFunc object itself is named 'this' in the Expression
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
%                     The Constraint may also be defined as a structure with
%                       min: a vector of minimal values per parameter
%                       max: a vector of maximal values per parameter
%                       fixed: a vector of 0(free) and 1(fixed) per parameter
%                       set: a cell of expressions per parameter
%                       Expression: an Expression to be evaluated
%                         The parameter names surrounded by "" are replaced
%                         by the corresponing p(n)
%
% Creating the object:
% ------------------------------------------------------------------------------
%   From a character string
%     the Expression should make use of x,y,z,t,u,v,w to denote axes of rank 1-7,
%     and the model Parameters are specified using 'p(n)' vector elements.
%       iFunc('p(1)*x+p(2)');
%
%   From a structure, with iFunc object fields (see above) and alias fields:
%     x0          -> Guess
%     objective   -> Expression
%     pars        -> Parameters
%     function    -> Name
%
%   From a function handle
%     The function should have syntax model(p, x,y,...) and return the model value.
%       iFunc(@(p,x) p(1)*x+p(2))
%
%   From a script/function 
%     The file should use 'p' as adjustable parameters, and optionally axes 'x', 'y', ...
%
%   From a SpinW object
%       iFunc(sw_model('squareAF',2,0))
%
%   From a file in MAT, M, YAML, JSON, XML format holding a model (iFunc/save)
%       iFunc('filename_model'). See iFunc/save
%
%   From a data set (iData) or data file
%       iFunc(iData('filename')); parameters allow to scale and shift.
%
% Using the object:
% ------------------------------------------------------------------------------
%   Once the object has been created, you can evaluate it with: object(p, x, y, ...)
%   The usual mathematical operators can be used to manipulate iFunc objects.
%   You can use guessed Model parameters with syntax: object('guess', x,y, ...)
%
%   The syntax iFunc(iData object) evaluates the iFunc model using the iData
%     object axes, and returns the model value as a numerical array. To get the
%     result as an iData, use syntax iData(iFunc object) [the opposite syntax].
%
%   The Model expression is               object.Expression
%   The Model parameter names are         object.Parameters
%   The Model parameter values (when set) object.ParameterValues or object.p
%   A single parameter can be accessed    object.<ParameterName> or object.p(index)
%
% Scan model parameters:
% ------------------------------------------------------------------------------
% It is possible to scan model parameters when using vectors as value for the
% parameters. To achieve that, the parameter values must be given as a named 
% structure.
% For instance, to scan the Intensity parameter in a gaussian, use:
%   model = gauss1;
%   p.Intensity= [ .5:.25:2 ];      % from 0.5 to 2 by steps of .25
%   v=iData(model, p);              % evaluate as iData sets.
%   plot(cat(2,v));                 % plot the scan as a surface
%
% Optimise model parameters:
% ------------------------------------------------------------------------------
% To optimise model parameters, you should first fix the non-varying
% parameters, and possibly bound the others. Then the optimisation is launched with
% any optimiser. To maximise the model, use '-model' as argument, as in the example:
%   model = gauss1;
%   fix(model, 'all'); model.Intensity='free';
%   model.Intensity=[0 1 2];    % bounds and starting value
%   model.HalfWidth = 0.5;
%   p = fmin(-model, [])        % return the optimal parameters to maximise the model
%
% Type <a href="matlab:doc(iFunc)">doc(iFunc)</a> to access the iFit/iFunc Documentation.
% iFunc is part of iFit http://ifit.mccode.org 
% 
% input:  s: iFunc, string, structure, function handle, spinw
% output: b: model object (iFunc)
% ex:     b=iFunc(@(p,x) p(1)*x+p(2)); 
%         b=iFunc('p(1)*x+p(2)');
%         b=iFunc('signal=p(1)*x+p(2);');
%         plot(x, b(p, x))
%         b.p=p; plot(b)
%
% Version: $Date$
% See also iFunc, iFunc/feval, iFunc/plot, iFunc/fits, iFunc/save, methods
% (c) E.Farhi, ILL. License: EUPL.

properties  % create the empty iFunc object structure
  % create a new iFunc object
  Tag         = 0;
  Date        = clock;  % creation date
  Name        = ''; % the function Name
  Description = ''; % the Description
  Parameters  = {}; % the function parameters, a cellstr or words or a single
                      % string of words separated by spaces, or a structure with parameter values
  Guess       = 'automatic'; % the default parameter set or expression or function handle
                      % char or function_handle(x,y,..., signal, Parameters{})
                      % or 'automatic'
  Expression  = ''; % the expression to evaluate to get the function value
                      % char or function_handle(p, x,y, ...)
  Constraint  = ''; % code to evaluate before computing the Expression
                      % can be: a function_handle
                      %         a char/cellstr
                      %         a vector (length p)
                      %         a structure with fields: min, max, fixed
  Dimension   = 0;  % function dimensionality (1,2,3,4...) 0=scalar=empty.
                      % a negative dimension is used to indicate a variable dimensionality.
  ParameterValues = [];
  Eval        = ''; % code to evaluate for the model value
  UserData    = '';
  Duration    = 0;
end % properties  

methods
  function a=iFunc(varargin)
  
    persistent id % the id number for all new objects
    
    % assign the ID and Tag
    if isempty(id),   id=0; end
    if id > 1e6, id=0; end
    if id <=0, 
      id = a.Date;
      id = fix(id(6)*1e4); 
    else 
      id=id+1;
    end

    a.Tag      = [ 'iF' sprintf('%0.f', id) ];

    if nargin == 0, return; end
    if length(varargin) > 1 % import data to create the object array
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
    elseif isa(varargin{1}, 'iData')
      % create a model from a fixed data set
      
      a = iFunc_from_iData(varargin{1});
      
    else   % import data to create a single object
      % can be a structure, or an expression, or a function_handle
      this = varargin{1};
      
      if ischar(this) && ~isempty(dir(this))  % a filename
        [p,f,e] = fileparts(this);
        if strcmp(lower(e), '.m')
          try
            run(this);
            if isa(ans, 'iFunc'), a=ans; end
          end
          if isempty(a)
            try
              a = eval(f);
            end
          end
        elseif strcmp(lower(e), '.mat')
          this=load(this);
          f=fieldnames(this);
          for index=1:numel(f)
            if isa(this.(f{index}), 'iFunc')
              a = [ a this.(f{index}) ];
            end
          end
        end
        if isempty(a) % other: JSON, YAML, XML will generate a structure when read
          a = iFunc(iLoad(this));
        end
      elseif ischar(this) % ----------------------------------------------------------
        % check if this is a predefined function
        isfunction = 0;
        try
          if exist(this) == 2 && abs(nargin(this)) > 0 && nargout(this) == 1
            isfunction = 1;
          end
        end
        if isfunction
          a=iFunc(feval(this));
          a.Name = this;
        else
          % from a string/expression: analyse the expression to get the
          % dimension and parameter names
          a=iFunc;
          if exist(this, 'file') || exist([ this '.m' ], 'file')
            a.Name = this;
            if exist([ this '.m' ], 'file'), this = [ this '.m' ]; end
            disp([ mfilename ': building model from file ' this ' content (script/function using p,x,y,...).']);
            this   = fileread(this);
          end
          a.Expression = this;
        end
      elseif iscell(this)
        for index=1:numel(this)
          try 
            a = [ a iFunc(this{index}) ];
          end
        end
        return
      elseif isa(this, 'iFunc') % ------------------------------------------------
        % make a copy of the initial object
        a = this;
      elseif isa(this, 'sw') | isa(this, 'spinw') % ------------------------------------------------
        % convert spinw object into iFunc
        a = sqw_spinw(this);
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
        if isfield(this, 'Dimension'),       a.Dimension=this.Dimension; end
        if isfield(this, 'ParameterValues'), a.ParameterValues=this.ParameterValues; end
        if isfield(this, 'UserData'),        a.UserData=this.UserData; end
        
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
      
    end % if nargin
  end % iFunc init
end % methods

end % classdef

% ------------------------------------------------------------------------------



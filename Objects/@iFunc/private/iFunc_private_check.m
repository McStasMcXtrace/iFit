function a = iFunc_private_check(a)
% iFunc_private_checkexpr: check the function expression/object

  if numel(a) > 1
    b = [];
    for index=1:numel(a)
      b = [ b iFunc_private_check(a(index)) ];
    end
    b = reshape(a, size(a));
    return
  end

  nb_pars             = 0;
  dim                 = 0;
  
  % Constraint can be of char, cellstr, function_handle, scalar, vector, structure
  const.min   = nan*ones(length(a.Parameters),1);
  const.max   = nan*ones(length(a.Parameters),1);
  const.fixed = zeros(length(a.Parameters),1);
  const.set   = cell(length(a.Parameters),1);
  const.eval  = '';
  if ~isstruct(a.Constraint)
    % build a structure of Constraints: these are evaluated in feval
    if ischar(a.Constraint) || iscellstr(a.Constraint) || isa(a.Constraint,'function_handle')
      const.eval = char(a.Constraint);          % iFunc.Constraint = char or cellstr or fhandle
    elseif isnumeric(a.Constraint)
      if length(a.Constraint) == length(a.Parameters) % iFunc.Constraint = vector
        const.fixed = a.Constraint(:);
      elseif length(a.Constraint) == 1                % iFunc.Constraint = scalar (0 or 1)
        const.fixed = const.fixed*a.Constraint;
      elseif ~isempty(a.Constraint)
        error(['iFunc:' mfilename ], [mfilename ': the model ' a.Tag ' Constraint should be scalar or vector of length ' ...
          num2str(length(a.Parameters)) ' (Parameters).' ]);
      end
    else
      error(['iFunc:' mfilename ], [mfilename ': the model ' a.Tag ' Constraint should be a char or cellstr, function_handle, scalar or vector, but not a ' ...
        class(a.Constraint) '.' ]);
    end
    a.Constraint = const;
  else
    % update the Constraint structure, search for min, max, fixed, Expression
    f=fieldnames(const);
    for index=1:length(f)
      if isfield(a.Constraint,f{index}) %  the Constraint given by user has the field from const
        v_new = a.Constraint.(f{index});
        l     = length(v_new);
        v_old = const.(f{index});
        if isnumeric(v_old) && l <= length(a.Parameters)  % this is a numeric field
          v_old(1:l) = v_new;
          const.(f{index}) = v_old; % update the values
        elseif ~isempty(v_new)
          if numel(v_new) < numel(v_old)
            v_new(numel(v_new)+1:numel(v_old)) = v_old(numel(v_new)+1:numel(v_old));
          end
          const.(f{index}) = v_new;
        end
      end
    end
    a.Constraint = const;
  end
  
  % Expression can be of char, cellstr and function_handle
  if isa(a.Expression,'function_handle')
    expr = func2str(a.Expression);
  elseif ischar(a.Expression) || iscellstr(a.Expression) 
    expr = char(a.Expression);
  else
    error(['iFunc:' mfilename ], [mfilename ': the model ' a.Tag ' Expression should be a char or function_handle or cellstr, not a class ' ...
      class(a.Expression) '.' ]);
  end

  % handle multiple lines in Expression -> single char with \n at end of lines (except last)
  n_expr = '';
  for index=1:size(expr, 1)
    d = strtrim(expr(index,:));
    if isempty(d), continue; end
    if d(end) ~= ';', d = [ d ';' ]; end
    if index == size(expr, 1)
      n_expr = [ n_expr d ];
    else
      n_expr = [ n_expr d sprintf('\n') ];
    end
  end
  expr = n_expr;
  
  if ischar(a.Parameters)
    pars              = strread(a.Parameters, '%s','delimiter',' ;,''"{}'); % is a cellstr
  elseif isstruct(a.Parameters)
    pars              = fieldnames(pars);
  elseif iscellstr(a.Parameters)
    pars              = a.Parameters;
  end
  
  if ~isempty(pars) && ~iscellstr(pars)
    warning([mfilename ': the model ' a.Tag ' parameters should be a char or structure or cellstr, not a class ' ...
      class(pars) '. Setting new parameter names.' ]);
    pars = {};
  end

  if strncmp(a.Guess, 'auto',4), a.Guess = ''; end
  
  % chech for function_handles stored as chars in Guess and Expression
  m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
  m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));
  if ischar(a.Guess) && numel(a.Guess) > 1
    if a.Guess(1) == '@', a.Guess = str2func(a.Guess); end
  end
  if ischar(a.Expression) && numel(a.Expression) > 1
    if a.Expression(1) == '@', a.Expression = str2func(a.Expression); end
  end

  if ~isempty(expr)
    % first count the initial p(n) fields and get the used indices
    % compute the number of parameters and corresponding parameter names
    % regexp: p(...) with '[0-9:]' characters inside
    nb_pars = [ regexp(expr,'\<p\(([\[0-9\:,;\s\]]+)\)','tokens')  ...
                regexp(expr,'_p\(([\[0-9\:,;\s\]]+)\)','tokens') ]; % return the tokens
    if ~isempty(nb_pars)
      % assemble all matches
      n = '';
      for index=1:length(nb_pars)
        this = nb_pars{index};
        n = [ n this{1}  ' ' ];
      end
      nb_pars  = unique(str2num([ '[' n ']' ]));
      used_pars= zeros(1,max(nb_pars));
      used_pars(nb_pars) = 1;
    else
      used_pars=[];
      nb_pars=[];
    end

    if ~isempty(nb_pars), 
      nb_pars = max(nb_pars);
    else nb_pars=0; end
    NL = sprintf('\n');
    e = textwrap(cellstr(char(a.Expression)),80);
    if length(e) > 6, e=[ e(1:3); '...'; e((end-2):end) ]; end
    e = sprintf('%s\n', e{:});
    if ~isempty(pars)
      % check if the number of parameters used in the expression matches the parameter names
      if nb_pars && length(pars) ~= nb_pars
        fprintf(1, [ '%s: Warning: iFit model ''%s''' NL ...
                     '    ''%s''' NL ...
                     '  The apparent number of parameters %d used in the expression p(1:%d) is not' NL ...
                     '  the same as the number of parameter names %d. Using the parameter names as valid input.\n' ], ...
                     mfilename, a.Tag, e, nb_pars, nb_pars, length(pars));

        % pars = []; % will guess parameter names
      end
      nb_pars = length(pars);
    end
    
    if ~isempty(a.Guess) && ~isa(a.Guess, 'function_handle')
      % check if the number of parameters used in the expression matches the parameter default values
      if ~ischar(a.Guess) && isnumeric(a.Guess) && numel(a.Guess) ~= nb_pars
        fprintf(1, [ '%s: Warning: iFit model ''%s''' NL ...
                     '    ''%s''' NL ...
                     '  The number of parameters %d used in the expression' NL ...
                     '  may not be the same as the number of default parameter values %d\n' ], ...
                     mfilename, a.Tag, e, nb_pars, numel(a.Guess));
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
              if length(regexp(expr, '\<v\>'))
                dim = 6;
                if length(regexp(expr, '\<w\>'))
                  dim = 7;
                end
              end
            end
          end
        end
      end
    end % if 'x' (1d)
  end % if ‪~empty(expr)

  % return when model can not be defined
  if isempty(expr) && ~dim
    warning('iFunc:emptyModel', [ 'model ' a.Tag ' Expression does not contain any axes (x,y,z...). Constant model ?' ]); 
  end
  if nb_pars == 0 && ~dim
    warning('iFunc:emptyModel', [ 'model ' a.Tag ' Expression does not contain any parameter (p). Constant model ?' ]); 
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
    a.Expression = cellstr(e);
  end % else keep function handle value

  if a.Dimension <= 0 && dim > 0, a.Dimension  = dim; end

  a.Parameters      = pars;
  a = iFunc_private_check_parnames(a, nb_pars);
  
  % check parameter values (if any defined yet): should be as many as parameter names
  val = [];
  if isstruct(a.Parameters) && length(fieldnames(a.Parameters)) == length(pars)
    val = cell2mat(struct2cell(a.Parameters));
  elseif ~isempty(a.ParameterValues)
    val = a.ParameterValues;
  end
  
  if ~isempty(val)
    if numel(pars) ~= numel(val)
      fprintf(1,'%s: model %s "%s": the number of model parameter values %d does not match the number of parameter names %d.\n', mfilename, a.Tag, a.Name, numel(val), length(pars));
      val = [ val(:) ; zeros(size(pars(:))) ];
      val = val(1:numel(pars));
    end
    a.ParameterValues = val;
  end
  
  a.Eval = cellstr(a);
  a.Eval = a.Eval(:);
  
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

% end iFunc_private_check

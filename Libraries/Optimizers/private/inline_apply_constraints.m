function pars = inline_apply_constraints(pars, constraints, options)
  % take into account constraints on parameters, and perform stop condition checks
  if isfield(constraints, 'step') % restrict parameter change
    parsStep    = pars(:) - constraints.parsPrevious(:);
    index       = find(constraints.steps(:) & abs(parsStep) > abs(constraints.steps(:)) & ~isnan(constraints.steps(:)));
    if ~isempty(index), 
      parsStep    = sign(parsStep).*abs(constraints.steps(:));
      pars(index) = constraints.parsPrevious(index) + parsStep(index);
    end
  end
  if isfield(constraints, 'min')    % lower bound for parameters
    index = find(pars(:) < constraints.min(:) & isfinite(constraints.min(:)));
    if ~isempty(index), pars(index) = constraints.min(index); end
  end
  if isfield(constraints, 'max')    % upper bound for parameters
    index = find(pars(:) > constraints.max(:) & isfinite(constraints.max(:)));
    if ~isempty(index), pars(index) = constraints.max(index); end
  end
  if isfield(constraints, 'fixed')  % fix some parameters
    index = find(constraints.fixed & ~isnan(constraints.fixed));
    if ~isempty(index), pars(index) = constraints.parsStart(index); end
  end
  if isfield(constraints, 'eval')   % evaluate expression with 'p'
    p = pars;
    try
      if isa(constraints.eval, 'function_handle')
        p = feval(constraints.eval, p);
      elseif ischar(constraints.eval)
        eval([ constraints.eval ';' ]);
      end
    end
    pars = p;
  end

  pars=pars(:)'; % parameters is a row
end % inline_apply_constraints

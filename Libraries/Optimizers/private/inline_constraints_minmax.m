function [constraints, constraints_var] = inline_constraints_minmax(pars, constraints)
% define default min max in constraints, needed by bounded optimizers
  if ~isfield(constraints, 'min')
    constraints.min = NaN*ones(size(pars));
  end
  for i=find(isnan(constraints.min(:)'));
    constraints.min(i) = -2*abs(pars(i)); % default min values
    if pars(i) == 0
      constraints.min(i) = -1;
    end
  end
  if ~isfield(constraints, 'max')
    constraints.max = NaN*ones(size(pars));
  end
  for i=find(isnan(constraints.max(:)'));
    constraints.max(i) =  2*abs(pars(i)); % default max values
    if pars(i) == 0
      constraints.max(i) = 1;
    end
  end
  % return as well only the variable parameters constraints
  constraints_var.min=constraints.min(constraints.index_variable);
  constraints_var.max=constraints.max(constraints.index_variable);
  if isfield(constraints,'step')
    constraints_var.step=constraints.step(constraints.index_variable);
  end
  if isfield(constraints,'set')
    constraints_var.set=constraints.set(constraints.index_variable);
  end
end % inline_constraints_minmax

function [c, pars] = inline_objective_all(fun, pars, options, constraints, varargin)
  % criteria to minimize, fun returns a scalar, or vector which is summed
  if nargin < 3, varargin={}; end
  % apply constraints on pars first
  pars                = inline_apply_constraints(pars,constraints,options); % private function
  % compute criteria
  c = feval(fun, pars, varargin{:});         % function=row vector, pars=column
  c = double(c(:)');
  switch options.optimizer
  case {'fminlm','LMFsolve'}
    % LMFsolve supports criteria as a vector of residuals, which sum is the criteria
    % but gradient is used to guide the optimizer
    if length(c) == 1,
      c = c*ones(1,10)/10;
    end
  otherwise
    c = sum(c);
  end
  
end % inline_objective_all

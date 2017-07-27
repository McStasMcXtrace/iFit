function result=test_Optimizers

  o = fits(iData); % get all Optimizers
  
  banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2; % Rosenbrock
  p = [];
  for index=1:length(o)
    [x,fval] = feval(o{index}, banana,[-1.2, 1]); % solution is x=[1 1]
    p = [ p ; x ];
  end
  
  % now test with a struct input
  [x1,fval] = feval('fminimfil', banana,struct('A',-1.2, 'B', 1));
  
  % now test with a constraint
  [x2,fval] = feval('fminimfil', banana,[-1.2 1],'', [0 1]);
  
  % now test with a constraint and named parameters
  [x3,fval] = feval('fminpso', banana,struct('A',-1.2, 'B', 1),'', [0 1]);

  if all(mean(abs(p)) > 0.5) ...
    && isstruct(x1) && abs(x1.A)-1 < .2 && abs(x1.B)-1< .2 ...
    && all(abs(x2)-1 < .2) ...
    && isstruct(x3) && abs(x3.A) -1 < .2
    result = [ 'OK     ' mfilename ' (' num2str(numel(o)) ' optimizers)' ];
  else
    result = [ 'FAILED ' mfilename ];
  end

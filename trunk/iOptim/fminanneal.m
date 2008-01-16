function [pars,fval,message,output] = fminanneal(fun, pars, options)
% [MINIMUM,FVAL,MESSAGE,OUTPUT] = FMINANNEAL(FUN,PARS,[OPTIONS]) simulated annealing optimizer
%
%  Simulated annealing optimizer. Prameters are guessed randomly with a decreasing
%  noise.
%
% Input:
%  FUN is a function handle (anonymous function or inline) with a loss
%  function, which may be of any type, and needn't be continuous. It does,
%  however, need to return a single value.
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the simulated annealing, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fminanneal')
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          MESSAGE return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
%  Contrib:
%    Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
%    Simulated Annealing. _Science, 220_, 671-680.
%   joachim.vandekerckhove@psy.kuleuven.be 2006/04/26 12:54:04
%
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-4;
  options.TolX   =1e-5;
  options.MaxIter=500;
  options.MaxFunEvals=10000;
  pars = options;
  return
end

if nargin <= 2
	options=[];
end
if isempty(options)
  options=feval(mfilename, 'defaults');
end

% call the optimizer
[pars,fval,message,output] = anneal(fun, pars, options);

% private function ------------------------------------------------------------

function [parent,fval,message,output] = anneal(loss, parent, options)

% main settings
%newsol      =@(x) (x+(randperm(length(x))==length(x))*randn/100);      % neighborhood space function
Tinit       = 1;        % initial temp
minT        = 1e-8;         % stopping temp
cool        = @(T) (.8*T);        % annealing schedule
minF        = options.TolFun;
max_consec_rejections = options.MaxFunEvals;
max_try     = options.MaxIter;
max_success = 20;
switch options.Display
case 'final'
  report=1;
case 'iter'
  report=2;
otherwise
  report=0;
end
k = 1;                           % boltzmann constant

% counters etc
iterations     = 0;
success  = 0;
consec   = 0;
T        = Tinit;
initenergy = feval(loss, parent);
fval     = initenergy;
funcount    = 0;
if report==2, fprintf(1,'\n  T = %7.5f, %s = %10.5f\n',T,func2str(loss),fval); end
message  = '';

while 1
  current = parent; 

  % % Stop / decrement T criteria
  if iterations >= max_try || success >= max_success;
      if consec >= max_consec_rejections
          message = 'Maximum consecutive rejections exceeded (options.MaxFunEvals)';
          break;
      elseif T < minT
          message = 'Minimum temperature reached';
          break;
      else
          T = cool(T);  % decrease T according to cooling schedule
          if report==2, % output
              fprintf(1,'Iter=%5i  T = %7.5f, %s = %10.5f\n',iterations,T,func2str(loss),fval);
          end
          iterations = iterations+1; % just an iteration counter
          success = 1;
      end
  end
%this inline is much more faster than original anneal.m version
% newparam = newsol(current);
  newparam = current;
  rand_ind = ceil(rand*length(current));
  stdx=0.05;
  newparam(rand_ind) = newparam(rand_ind)+randn*stdx;
  funcount = funcount + 1;
  newfval = feval(loss,newparam);
  
  % std stopping conditions
  [istop, message] = fmin_private_std_check(newparam, newfval, iterations, funcount, options);
  
  if istop
    parent = newparam; 
    fval = newfval;
    break
  end

  if (fval-newfval > min(options.TolFun/2, 1e-6))
      parent = newparam;
      fval = newfval;
      success = success+1;
      consec = 0;
  else
      if (rand < exp( (fval-newfval)/(k*T) ));
          parent = newparam;
          fval = newfval;
          success = success+1;
      else
          consec = consec+1;
      end
  end

end % while

if report;
    fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
    fprintf(1, '  Final temperature:       \t%g\n', T);
    fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
    fprintf(1, '  Number of function calls:\t%i\n', funcount);
    fprintf(1, '  Minimal final %s:        \t%g\n', func2str(loss), fval);
end

output.iterations = iterations;
output.algorithm  = mfilename;
output.message    = message;
output.funcCount  = funcount;


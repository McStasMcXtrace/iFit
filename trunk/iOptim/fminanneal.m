function [pars,fval,exitflag,output] = fminanneal(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINANNEAL(FUN,PARS,[OPTIONS]) simulated annealing optimizer
%
%  Simulated annealing minimization. Parameters are scanned randomly with a decreasing
%  noise.
% 
% Calling:
%   fminanneal(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminanneal(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminanneal(banana,[-1.2, 1])
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
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Reference:
%    Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
%    Simulated Annealing. _Science, 220_, 671-680.
% Contrib:
%   joachim.vandekerckhove@psy.kuleuven.be 2006/04/26 12:54:04
%
% Version: $Revision: 1.5 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-4;
  options.TolX   =0;
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

options=fmin_private_std_check(options);

options.algorithm  = [ 'Simulated Annealing (be Vandekerckhove) [' mfilename ']' ];

% call the optimizer
[pars,fval,exitflag,output] = anneal(fun, pars, options);

% private function ------------------------------------------------------------

function [parent,fval,istop,output] = anneal(loss, parent, options)

% main settings
%newsol      =@(x) (x+(randperm(length(x))==length(x))*randn/100);      % neighborhood space function
Tinit       = 1;        % initial temp
minT        = 1e-8;         % stopping temp
cool        = @(T) (.8*T);        % annealing schedule
minF        = options.TolFun;
max_consec_rejections = options.MaxFunEvals;
max_try     = options.MaxIter;
max_success = 20;
k = 1;                           % boltzmann constant

% counters etc
iterations     = 0;
success  = 0;

consec   = 0;
T        = Tinit;
initenergy = feval(loss, parent);
fval     = initenergy;
funcount    = 0;
istop    =0;

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, loss, parent, fval);
end

while 1
  current = parent; 

  % % Stop / decrement T criteria
  if iterations >= max_try || success >= max_success;
      if consec >= max_consec_rejections
          message = 'Maximum consecutive rejections exceeded (options.MaxFunEvals)';
          istop=-7;
          break;
      elseif T < minT
          message = 'Minimum temperature reached';
          istop=-8;
          break;
      else
          T = cool(T);  % decrease T according to cooling schedule
          iterations = iterations+1; % just an iteration counter
          success = 1;
          if strcmp(options.Display, 'iter')
            fmin_private_disp_iter(iterations, funcount, loss, newparam, newfval);
          end
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
  [istop, message] = fmin_private_std_check(newparam, newfval, iterations, funcount, options, parent, fval);
  
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

% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = iterations;
output.algorithm  = options.algorithm;
output.message    = message;
output.funcCount  = funcount;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, loss, parent, fval);
end



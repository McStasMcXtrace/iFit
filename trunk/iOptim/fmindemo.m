function fmindemo(dim, verbose)
% fmindemo(dimensionality=vector, verbose=0|1)
%
% ex: fmindemo([2 10]);

if nargin == 0
  fmindemo(2);
  fmindemo(10);
  fmindemo(100);
  return
end

if nargin == 1
  verbose=0;
end

if isnumeric(dim) & length(dim) > 1
  for i=1:length(dim)
    fmindemo(dim(i), verbose);
  end
  return
end

optimizers={...
  'fminanneal', 'fminrand', ...
  'fmincmaes','fminga', ...
  'fmingradrand','fminpowell','fminralg', ...
  'fminkalman',...
  'fminsimplex','fminsearchOS',...
  'fminswarm','fminswarmhybrid'};

problems={...       % Unimodal functions
  'fsphere', ...
  'fsphereoneax', ...
  'frandsphere', ...
  'fspherehull', ...
  'fneumaier3', ...
  'flogsphere', ...
  'fexpsphere', ...
  'fbaluja', ...
  'fschwefel', ...
  'fplane', ...
  'ftwoaxes', ...
  'fparabR', ...
  'fsharpR', ...
  'frosen', ...
  'fschwefelrosen1', ...
  'fschwefelrosen2', ...
  'fackley', ...    % Multimodal functions
  'fconcentric', ...
  'fgriewank', ...
  'fschaffer', ...
  'fschwefelmult', ...
  'ftwomaxtwo'};

dimensionality=randn(1,dim);

w=warning;
warning off


startpars = dimensionality;
for index_problem=1:length(problems)
  fprintf(1, 'Problem: %15s [%i parameters] ', ...
    problems{index_problem}, length(dimensionality) );
  %    disp(startpars)
  fun=str2func(problems{index_problem});
  for index=1:length(optimizers)
    fprintf(1,'%i', index);
    t0 = clock;
    options=feval(optimizers{index},'defaults');
    options.Display='off';
    maxit = options.MaxIter; if ischar(maxit), maxit=str2num(maxit); end
    options.MaxIter=max(500, maxit);
    maxfn = options.MaxFunEvals; if ischar(maxfn), maxfn=str2num(maxfn); end
    options.MaxFunEvals=max(5000, maxfn);
    try
      [pars, fval, flag, out] = feval(optimizers{index}, fun, startpars(:), options);
    catch
      pars=[]; fval=Inf; out.funCount=Inf; flag=-4;
      out.iterations=Inf; out.algorithm=[ 'ERROR: ' optimizers{index} ];
      lasterr
    end
    duration(index)=etime(clock,t0);
    funcount(index)=out.funcCount;
    iterate(index) =out.iterations;
    if iterate(index)>=options.MaxIter, f='nITR';
    elseif funcount(index)>=options.MaxFunEvals, f='nFUN';
    elseif fval <= options.TolFun, f='CFUN';
    else
      switch flag
      case 0, f='COK ';
      case -1,f='CFUN';
      case -2,f='nITR';
      case -3,f='nFUN';
      case -4,f='*INF';
      case -5,f='CX  ';
      otherwise
        if flag<0, f='n'; else f='C'; end
      end
    end
    fvals(index)   =fval;
    flags{index}   =f;
    algo{index}    =out.algorithm;
    abstract{index_problem, index} = f(1);
    fprintf(1,'%s ', f(1));
  end
  disp(' ');
  % sort results
  [dummy, sorti] = sort(duration.*fvals.*fvals);
  
  % compute median duration
  mtime = mean(duration)+std(duration);
  for index=1:length(optimizers)
    f = abstract{index_problem, index};
    if f(1) == 'C',   
      if duration(index) > mtime, f='-'; 
      else f='+'; end
    elseif f(1)=='n', f=' '; 
    else f='*'; end
    abstract{index_problem, index} = f(1);
  end

  if verbose
    % display results
    fprintf(1, 'Index    Time FuncEval Iter.      Fval Status Algo\n');
    for j=1:length(optimizers)
      index=sorti(j);
      fprintf(1, '%5i %7.3g %7i %5i %10.2g %4s   %s\n', ...
        index, duration(index), funcount(index), iterate(index), fvals(index), flags{index}, algo{index});
    end
  end
end

% display for given dimensionality an abstract of tests
% matrix: rows=tests columns=optimizers, values=sorting index
fprintf(1, '               ');
for j=1:length(optimizers)
  opt=optimizers{j};
  opt=opt(5:end);
  fprintf(1, '%3s ', opt(1:min(3, length(opt)))); 
end
fprintf(1, '\n');
for i=1:length(problems)
  fprintf(1, '%15s ', problems{i});
  for j=1:length(optimizers)
    fprintf(1, '%s   ', abstract{i, j});
  end
  fprintf(1, '\n');
end

warning(w);



% TEST functions taken from CMA-ES

function f=fsphere(x)

  f=sum(x.^2);

function f = fsphereoneax(x)
  f = x(1)^2;
  f = mean(x)^2;
  f=abs(f);

function f=frandsphere(x)
  N = size(x,1);
  idx = ceil(N*rand(7,1));
  f=sum(x(idx).^2);
  f=abs(f);

function f=fspherelb0(x, M) % lbound at zero for 1:M needed
  if nargin < 2 M = 0; end
  N = size(x,1);
  % M active bounds, f_i = 1 for x = 0
  f = -M + sum((x(1:M) + 1).^2);
  f = f + sum(x(M+1:N).^2);
  f=abs(f);

function f=fspherehull(x)
  % Patton, Dexter, Goodman, Punch
  % in -500..500
  % spherical ridge through zeros(N,1)
  % worst case start point seems x = 2*100*sqrt(N)
  % and small step size
  N = size(x,1);
  f = norm(x) + (norm(x-100*sqrt(N)) - 100*N)^2;
  f=abs(f);

function f=fellilb0(x, idxM, scal) % lbound at zero for 1:M needed
  N = size(x,1);
  if nargin < 3 || isempty(scal)
    scal = 100;
  end
  scale=scal.^((0:N-1)/(N-1));
  if nargin < 2 || isempty(idxM)
    idxM = 1:N;
  end
  %scale(N) = 1e0;
  % M active bounds
  xopt = 0.1;
  x(idxM) = x(idxM) + xopt;
  f = scale.^2*x.^2;
  f = f - sum((xopt*scale(idxM)).^2);
%  f = exp(f) - 1;
%  f = log10(f+1e-19) + 19;

  f = f + 1e-19;
  f=abs(f);

function f=fsectorsphere(x, scal)
%
% This is deceptive for cumulative sigma control in large dimension:
% The strategy (initially) diverges for N=50 and popsize = 150.  (Even
% for cs==1 this can be observed for larger settings of N and
% popsize.) The reason is obvious from the function topology.
% Divergence can be avoided by setting boundaries or adding a
% penalty for large ||x||. Then, convergence can be observed again.
% Conclusion: for popsize>N cumulative sigma control is not completely
% reasonable, but I do not know better alternatives.
%
  if nargin < 2 || isempty (scal)
    scal = 1e3;
  end
  f=sum(x.^2);
  idx = find(x<0);
  f = f + (scal-1)^2 * sum(x(idx).^2);
  f=abs(f);

function f=fstepsphere(x, scal)
  if nargin < 2 || isempty (scal)
    scal = 1e0;
  end
  N = size(x,1);
  f=1e-11+sum(scal.^((0:N-1)/(N-1))*floor(x+0.5).^2);
  f=1e-11+sum(floor(scal.^((0:N-1)/(N-1))'.*x+0.5).^2);
  f=abs(f);
%  f=1e-11+sum(floor(x+0.5).^2);

function f=fstep(x)
  % in -5.12..5.12 (bounded)
  N = size(x,1);
  f=1e-11+6*N+sum(floor(x));
  f=abs(f);

function f=fneumaier3(x)
  % in -n^2..n^2
  % x^*-i = i(n+1-i)
  N = size(x,1);
%  f = N*(N+4)*(N-1)/6 + sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  f = sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  f=abs(f);

function f= flogsphere(x)
 f = 1-exp(-sum(x.^2));
 f=abs(f);

function f= fexpsphere(x)
 f = exp(sum(x.^2)) - 1;
 f=abs(f);

function f=fbaluja(x)
  % in [-0.16 0.16]
  y = x(1);
  for i = 2:length(x)
    y(i) = x(i) + y(i-1);
  end
  f = 1e5 - 1/(1e-5 + sum(abs(y)));
  f=abs(f);

function f=fschwefel(x)
  f = 0;
  for i = 1:size(x,1),
    f = f+sum(x(1:i))^2;
  end
  f=abs(f);

function f=fcigar(x, ar)
  if nargin < 2 || isempty(ar)
    ar = 1e3;
  end
  f = x(1)^2 + ar^2*sum(x(2:end).^2);
  f=abs(f);

function f=felli(x, lgscal, expon, expon2)
  % lgscal: log10(axis ratio)
  % expon: x_i^expon, sphere==2
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  if nargin < 2 || isempty(lgscal), lgscal = 3; end
  if nargin < 3 || isempty(expon), expon = 2; end
  if nargin < 4 || isempty(expon2), expon2 = 1; end

  f=((10^(lgscal*expon)).^((0:N-1)/(N-1)) * abs(x).^expon).^(1/expon2);
%  if rand(1,1) > 0.015
%    f = NaN;
%  end
f=abs(f);


function f=fellii(x, scal)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  if nargin < 2
    scal = 1;
  end
  f= (scal*(1:N)).^2 * (x).^2;
  f=abs(f);

function f=fplane(x)
  f=x(1);
  f=abs(f);

function f=ftwoaxes(x)
  f = sum(x(1:floor(end/2)).^2) + 1e6*sum(x(floor(1+end/2):end).^2);
  f=abs(f);

function f=fparabR(x)
  f = -x(1) + 100*sum(x(2:end).^2);
  f=abs(f);

function f=fsharpR(x)
  f = abs(-x(1)) + 30*norm(x(2:end));
  f=abs(f);

function f=frosen(x)
  if length(x) < 2 error('dimension must be greater one'); end
  f = 1e2*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);
  f=abs(f);

function f=frosenmodif(x)
  f = 74 + 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 ...
      - 400*exp(-sum((x+1).^2)/2/0.05);
  f=abs(f);

function f=fschwefelrosen1(x)
  % in [-10 10]
  f=sum((x.^2-x(1)).^2 + (x-1).^2);
  f=abs(f);

function f=fschwefelrosen2(x)
  % in [-10 10]
  f=sum((x(2:end).^2-x(1)).^2 + (x(2:end)-1).^2);
  f=abs(f);

%%% Multimodal functions

function f=fackley(x)
  % -32.768..32.768
  % Adding a penalty outside the interval is recommended,
  % because for large step sizes, fackley imposes like frand
  %
  N = size(x,1);
  f = 20-20*exp(-0.2*sqrt(sum(x.^2)/N));
  f = f + (exp(1) - exp(sum(cos(2*pi*x))/N));
  % add penalty outside the search interval
  f = f + sum((x(x>32.768)-32.768).^2) + sum((x(x<-32.768)+32.768).^2);
  f=abs(f);

function f = fbohachevsky(x)
 % -15..15
  f = sum(x(1:end-1).^2 + 2 * x(2:end).^2 - 0.3 * cos(3*pi*x(1:end-1)) ...
	  - 0.4 * cos(4*pi*x(2:end)) + 0.7);
	f=abs(f);

function f=fconcentric(x)
  % in  +-600
  s = sum(x.^2);
  f = s^0.25 * (sin(50*s^0.1)^2 + 1);
  f=abs(f);

function f=fgriewank(x)
  % in [-600 600]
  N = size(x,1);
  f = 1 - prod(cos(x'./sqrt(1:N))) + sum(x.^2)/4e3;
  % f = f + 1e4*sum(x(abs(x)>5).^2);
  % if sum(x(abs(x)>5).^2) > 0
  %   f = 1e4 * sum(x(abs(x)>5).^2) + 1e8 * sum(x(x>5)).^2;
  % end
  f=abs(f);

function f=frastrigin(x, scal, skewfac, skewstart, amplitude)
% by default multi-modal about between -30 and 30
  if nargin < 5 || isempty(amplitude)
    amplitude = 10;
  end
  if nargin < 4 || isempty(skewstart)
    skewstart = 0;
  end
  if nargin < 3 || isempty(skewfac)
    skewfac = 1;
  end
  if nargin < 2 || isempty(scal)
    scal = 1;
  end
  N = size(x,1);
  scale = 1;
  if N > 1
    scale=scal.^((0:N-1)'/(N-1));
  end
  % simple version:
  % f = amplitude*(N - sum(cos(2*pi*(scale.*x)))) + sum((scale.*x).^2);

  % skew version:
  y = scale.*x;
  idx = find(x > skewstart);
  if ~isempty(idx)
    y(idx) =  skewfac*x(idx);
  end
  f = amplitude * (N-sum(cos(2*pi*y))) + sum(y.^2);
  f=abs(f);

function f = fschaffer(x)
 % -100..100
  N = size(x,1);
  s = x(1:N-1).^2 + x(2:N).^2;
  f = sum(s.^0.25 .* (sin(50*s.^0.1).^2+1));
  f=abs(f);

function f=fschwefelmult(x)
  % -500..500
  N = size(x,1);
  %   f = - sum(x.*sin(sqrt(abs(x))));
  f = 418.9829*N - 1.27275661e-5*N - sum(x.*sin(sqrt(abs(x))));
  % penalty term
  f = f + sum(x(abs(x)>500).^2);
  f=abs(f);

function f=ftwomax(x)
  % Boundaries at +/-5
  N = size(x,1);
  f = -abs(sum(x)) + 5*N;
  f=abs(f);

function f=ftwomaxtwo(x)
  % Boundaries at +/-10
  N = size(x,1);
  f = abs(sum(x));
  if f > 30
    f = f - 30;
  end
  f = -f;
  f=abs(f);


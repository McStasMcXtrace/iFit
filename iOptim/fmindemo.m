function abstract=fmindemo(dim, verbose)
% fmindemo(dimensionality=vector, verbose=0|1)
%
% ex: fmindemo([2 10]);

if nargin == 0
  abstract{1} = fmindemo(2);
  abstract{2} = fmindemo(10);
  abstract{3} = fmindemo(100);
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
% ann bfg new imf cma  ga gra pow ral kal sim sea hoo swa swa ran mul 
optimizers={ ...
    'fminpowell','fminralg','fmincmaes','fminsimplex','fminnewton',...
    'fmingradrand', 'fminimfil', 'fminbfgs', ...
    'fminhooke', 'fminanneal', ...
    'fminga', ...
    'fminkalman',...
    'fminsearchOS', ...
    'fminswarm','fminswarmhybrid'};

if dim <=10
  optimizers{end+1} = 'fminrand';  % this one is too slow
  optimizers{end+1} = 'fminmulti'; % this one is not accurate
end

problems={...       % Unimodal functions
  'fsphere',      1, ...
  'fsphereoneax', 1, ...
  'frandsphere',  1, ...
  'fspherehull',  500, ...
  'fneumaier3',   5, ...
  'flogsphere',   1, ...
  'fexpsphere',   1, ...
  'fbaluja',      0.16, ...
  'fschwefel',    1, ...
  'fplane',       1, ...
  'ftwoaxes',     1, ...
  'fparabR',      1, ...
  'fsharpR',      1, ...
  'frosen',       1, ...
  'fschwefelrosen1', 10, ...
  'fschwefelrosen2', 10, ...
  'fackley', 32000, ...    % Multimodal functions
  'fconcentric', 600, ...
  'fgriewank', 600, ...
  'fschaffer', 100, ...
  'fschwefelmult', 500, ...
  'ftwomaxtwo', 10};

dimensionality=randn(1,dim);

w=warning;
warning off

disp([ 'Starting optimization demo with ' num2str(dim) ' parameters to find.' ]);
for index=1:length(optimizers)
  options=feval(optimizers{index},'defaults');
  opt=optimizers{index};
  if isfield(options, 'algorithm')
    alg=options.algorithm;
  else
    alg=opt;
  end
  if length(opt) > 5
    opt=opt(5:end);
    opt=opt(1:min(3, length(opt)));
  end
  fprintf(1, '%5s %s\n', opt, alg);
end

t1 = clock;

for index_problem=1:2:length(problems)
  fprintf(1, 'Problem: %15s [%i parameters] ', ...
    problems{index_problem}, length(dimensionality) );
  %    disp(startpars)
  startpars = dimensionality*problems{index_problem+1};
  fun=str2func(problems{index_problem});
  for index=1:length(optimizers)
    if ~verbose, fprintf(1,'%i', index); end
    t0 = clock;
    options=feval(optimizers{index},'defaults');
    %if verbose, options.Display='iter'; end
    maxit = options.MaxIter; if ischar(maxit), maxit=eval(maxit,'0'); end
    if isinf(maxit), maxit=0; end
    options.MaxIter=min(2000, max(250*dim, maxit));
    maxfn = options.MaxFunEvals; if ischar(maxfn), maxfn=eval(maxfn,'0'); end
    if isinf(maxfn), maxfn=0; end
    options.MaxFunEvals=min(20000, max(2500*dim, maxfn));
    options.TolX=0;
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
    if isnan(fval) | isinf(fval), fval=Inf; flag=-4; end
    if ~isinf(fval) & iterate(index)>=options.MaxIter, f='nITR';
    elseif ~isinf(fval) & funcount(index)>=options.MaxFunEvals, f='nFUN';
    elseif fval <= options.TolFun, flag=-1; f='CFUN';
    else
      switch flag
      case -1,f='Cfun';
      case -2,f='nitr';
      case -3,f='nFUN';
      case -4,f='*INF';
      case -5,f='Cx  ';
      case {-12, -9, -8}, f='.';
      otherwise
        if fval < options.TolFun, f='C';
        else f='n'; end
      end
    end
    fvals(index)   =fval;
    flags{index}   =f;
    algo{index}    =out.algorithm;
    params{index}=pars;
    abstract{index_problem, index} = f(1);
    if ~verbose, fprintf(1,'%s ', f(1)); end
  end
  disp(' ');
  % sort results
  criteria=duration.*fvals.*fvals;
  index=strmatch('C', flags);
  criteria(index) = criteria(index)*100; % put converged methods first
  [dummy, sorti] = sort(criteria);
  
  % compute median duration
  mtime = mean(duration)+std(duration);
  for index=1:length(optimizers)
    f = abstract{index_problem, index};
    if f(1) == 'C',   
      if duration(index) > mtime, f='-'; 
      else f='+'; end
    elseif f(1)=='.', f='.';
    elseif f(1)=='n', f=' '; 
    else f='*'; end
    abstract{index_problem, index} = f(1);
  end

  if verbose
    % display results
    fprintf(1, 'Index    Time FuncEval Iter.      Fval Status Pars      Algo\n');
    for j=1:length(optimizers)
      index=sorti(j);
      pars=params{index};
      if length(pars) > 5, pars=pars(1:5); end
      spars=mat2str(pars(:)');  % as a row
      if length(spars) > 50, spars=[ spars(1:47) ' ...' ]; end
      opt=optimizers{index};
      opt=opt(5:end);
      fprintf(1, '%5s %7.3g %7i %5i %10.2g %4s   %50s %s\n', ...
        opt(1:min(3, length(opt))), duration(index), funcount(index), iterate(index), fvals(index), flags{index}, spars, algo{index});
    end
  end
end

% display for given dimensionality an abstract of tests
% matrix: rows=tests columns=optimizers, values=sorting index
disp('Test results: +=fast; -=slower; .=local minimum; *=error; " "=does not converge');
fprintf(1, 't=%10.2g s ', etime(clock,t1));
for j=1:length(optimizers)
  opt=optimizers{j};
  opt=opt(5:end);
  fprintf(1, '%3s ', opt(1:min(3, length(opt)))); 
end
fprintf(1, '\n');
for i=1:2:length(problems)
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

function f=fspherehull(x)
  % Patton, Dexter, Goodman, Punch
  % in -500..500
  % spherical ridge through zeros(N,1)
  % worst case start point seems x = 2*100*sqrt(N)
  % and small step size
  N = size(x,1);
  f = norm2(x) + (norm2(x-100*sqrt(N)) - 100*N)^2;
  f=abs(f);

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
  f = abs(-x(1)) + 30*norm2(x(2:end));
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

function f=ftwomaxtwo(x)
  % Boundaries at +/-10
  N = size(x,1);
  f = abs(sum(x));
  if f > 30
    f = f - 30;
  end
  f = -f;
  f=abs(f);
  
function n=norm2(x)
x = x(:);
n=sqrt(sum(abs(x).*abs(x)));



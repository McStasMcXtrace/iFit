function [dp,covp,corp,jac,hessian] = inline_estimate_uncertainty(fun, pars, pars_all, options, constraints, varargin) %KP modified
% [dp,covp,corp] = inline_estimate_uncertainty(fun, pars, options)
% 
% Estimates the uncertainty around an optimization solution using
% the error matrix from the criteria jacobian inversion.
% 
% Calling:
%   [p,c,e,o]  = fmin(fun, p0, ...);
%   dp = o.parsHessianUncertainty;
%
% Input:
%  FUN is a function handle (anonymous function or inline) with a loss
%  function, which may be of any type, and needn't be continuous. It does,
%  however, need to return a single value.
%
%  PARS is a vector with initial guess of fitted parameters. You must input an
%  initial guess.
%
%  PARS_ALL is a vector with inicial guess of all parameters.
%
%  OPTIONS is a structure with settings for the simulated annealing, 
%  compliant with optimset. Default options may be obtained with
%     o=fmin_private_wrapper(optimizer,'defaults')
%
%  CONSTRAINTS is a vector specifying fixed parameters.
%
% Output:
%  DP is the gaussian uncertainty around PARS
%  COVP is the error matrix 
%  CORP is the correlation matrix
%  JAC  is the Jacobian
%  HESSIAN is the Hessian

  %DOSTAT NA VSTUP OBOJE - PARS A PARS_ALL, TAKY CONSTRAINTS

  %vector of indices of fitted parameters %KP modified
  ind = find(constraints.fixed == 0);
  
  n = length(pars);
  if nargin < 3, options=[]; end
  if nargin < 4, args={}; end
  % if isfield(options,'TolX') TolX = options.TolX; 
  % else TolX = 0; end

  % initialize the curvature matrix alpha = '1/2 d2 Chi2/dpi/dpj' (half Hessian)
  alpha= zeros(n);
  dp   = zeros(size(pars));
  val=feval(fun, pars_all, varargin{:}); %val=feval(fun, pars, varargin{:}); %KP modified
  mm=length(val); %number of measurement points 
  chisq= sum(val);
  
  covp = [];
  corp = [];
  jac  = [];
  hessian = [];
  %if TolX <= 0, 
  %  TolX = 0.01*pars;
  %end
  %if length(TolX) == 1
  %  dp   = TolX*ones(size(pars));
  %else
  %  dp   = TolX;
  %end
  %dp(find(dp == 0)) = 1e-5;
  dp=ones(size(pars_all))*1e-6;

  % we now build the error matrix 'alpha' and the Jacobian
  jac = zeros(n,1);
  
  for i=1:n %KP modified
    ii = ind(i);
    p = pars_all; p(ii) = p(ii)+dp(ii); chi1 = sum(feval(fun, p, varargin{:}));
    p = pars_all; p(ii) = p(ii)-dp(ii); chi2 = sum(feval(fun, p, varargin{:}));
    %p = pars; p(i) = p(i)+dp(i); chi1 = sum(feval(fun, p, varargin{:}));
    %p = pars; p(i) = p(i)-dp(i); chi2 = sum(feval(fun, p, varargin{:}));
    alpha(i,i) = (chi1-2*chisq+chi2)/2/dp(ii)/dp(ii); % diagonal terms
    jac(i) = (chi1-chisq)/dp(ii);
    
    for j=i+1:n
      jj = ind(j);
      p = pars_all; p(ii) = p(ii)+dp(ii); p(jj) = p(jj)+dp(jj); chi1 = sum(feval(fun,p, varargin{:}));
      p=pars_all; p(ii)=p(ii)+dp(ii); p(jj)=p(jj)-dp(jj); chi2=sum(feval(fun,p, varargin{:}));
      p=pars_all; p(ii)=p(ii)-dp(ii); p(jj)=p(jj)+dp(jj); chi3=sum(feval(fun,p, varargin{:}));
      p=pars_all; p(ii)=p(ii)-dp(ii); p(jj)=p(jj)-dp(jj); chi4=sum(feval(fun,p, varargin{:}));
      %p=pars; p(i)=p(i)+dp(i); p(j)=p(j)+dp(j); chi1=sum(feval(fun,p, varargin{:}));
      %p=pars; p(i)=p(i)+dp(i); p(j)=p(j)-dp(j); chi2=sum(feval(fun,p, varargin{:}));
      %p=pars; p(i)=p(i)-dp(i); p(j)=p(j)+dp(j); chi3=sum(feval(fun,p, varargin{:}));
      %p=pars; p(i)=p(i)-dp(i); p(j)=p(j)-dp(j); chi4=sum(feval(fun,p, varargin{:}));
      alpha(i,j)=(chi1-chi2-chi3+chi4)/8/dp(i)/dp(j);
      alpha(j,i)=alpha(i,j); % off diagonal terms (symmetric)
    end
  end
  if any(isnan(alpha(:))), return; end 
  hessian=2*alpha;
  alpha = alpha/chisq*(mm-n);      % normalized error matrix
  covp  = pinv(alpha);       % COV MATRIX
  dp    = sqrt(abs(diag(covp))); % uncertainty on parameters
  corp  = covp./(dp*dp');   % correlation matrix
  
end % inline_estimate_uncertainty

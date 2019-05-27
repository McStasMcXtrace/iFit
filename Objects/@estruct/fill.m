function c=fill(a, n)
% FILL fill-in over missing data.
%   Y = FILL(X) replaces the missing data in X by extra/interpolating
%   the non-missing elements. The non finite values (NaN or Inf) in X are
%   considered as missing data. The method is using a N-D discrete cosine
%   transform.
%   FILL uses an iterative process that converges toward the solution.
%
%   Y = FILL(X,N) uses N iterations. By default, N = 100. If you
%   estimate that FILL did not totally converge, increase N:
%   Y = FILL(X,1000);
%
%   For a faster fill-in, you may use RESIZE instead.
%
% Example: a=estruct([ ifitpath 'Data/Monitor_GV*']); b=hist(a); c=fill(b);
% Version: $Date$ $Version$ $Author$
% See also estruct, accumarray, hist, histc, estruct/plot, sum, estruct/interp, estruct/event, estruct/squeeze, estruct/pack

% private: inpaintn from Garcia D. 2012-2014

c=[];

if nargin < 2, n=[]; end
if isempty(n), n=100; end

% handle handle array as input
if numel(a) > 1
  c = zeros(estruct, numel(a), 1);
  for index=1:numel(a)
    c(index) = feval(mfilename, a(index), n);
  end
  return
end

% fill values=NaN
s      = get(a, 'Signal');
e      = subsref(a,struct('type','.','subs','Error'));
m      = subsref(a,struct('type','.','subs','Monitor'));
index0 = isnan(s);
if numel(m) == numel(s)
  index0= isnan(s) | m==0 | isnan(m);
end

if ~isempty(index0)
  s(index0) = NaN;
  s = inpaintn(s, n);
  if ~isscalar(e)
    e(index0) = NaN;
    e = inpaintn(e, n);
  end
  if ~isscalar(m)
    m(index0) = 1;
  end
else
  c = a; 
  return
end

% assemble final new object
c = zeros(a);

c = set(c, 'Signal',  s);
c = set(c, 'Error',   e);
c = set(c, 'Monitor', m);

history(c, mfilename, a, n);

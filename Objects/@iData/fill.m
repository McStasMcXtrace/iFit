function c=fill(a, n)
% FILL fill-in over missing data
%   Y = FILL(X) replaces the missing data in X by extra/interpolating
%   the non-missing elements. The non finite values (NaN or Inf) in X are
%   considered as missing data.
%
%   FILL uses an iterative process that converges toward the solution.
%   Y = FILL(X,N) uses N iterations. By default, N = 100. If you
%   estimate that FILL did not totally converge, increase N:
%   Y = FILL(X,1000);
%
% input:  a: object or array (iData)
%         n: number of iterations (integer)
% output: c: filled object (iData)
% ex:     a=iData([ ifitpath 'Data/Monitor_GV*']); b=hist(a); c=fill(b);
%
% Version: $Date$
% See also iData, accumarray, hist, histc, iData/plot, sum, iData/interp, iData/event, iData/squeeze, iData/pack

% private: inpaintn from Garcia D. 2012-2014

c=[];

if nargin < 2, n=[]; end
if isempty(n), n=100; end

% handle handle array as input
if numel(a) > 1
  c = zeros(iData, numel(a), 1);
  parfor index=1:numel(a)
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
  e(index0) = NaN;
  m(index0) = 1;
  s = inpaintn(s, n);
  e = inpaintn(e, n);
end

% assemble final new object
c = copyobj(a);

c = setalias(c, 'Signal',  s);
c = setalias(c, 'Error',   e);
c = setalias(c, 'Monitor', m);

function v = linspace(a,b,n)
% LINSPACE Linearly spaced objects.
%   LINSPACE(X1, X2) generates a row vector of 10 linearly
%   equally spaced objects between X1 and X2.
%   This corresponds to a 'morphing' from X1 to X2 in N steps.
%
%   LINSPACE(X1,X2,N) generates N objects between X1 and X2.
%
% Example: b=estruct(peaks); a=-b; c=linspace(a,b); numel(c) == 10
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/max, estruct/min, estruct/colon, estruct/logspace

if nargin <= 2
  n = [];
end
if nargin == 1, b=a; end

if isempty(n) || n <=0
  n=10;
end

if isa(a, 'estruct') && numel(a) > 1, a=a(1); end
if isa(b, 'estruct') && numel(b) > 1, b=b(end); end

% get scalar value (if specified as number of 1x1 object)
if     isempty(a), a=0; 
elseif isempty(b), b=0; 
elseif isscalar(a) & isa(a, 'estruct'), a=get(a,'Signal');
elseif isscalar(b) & isa(b, 'estruct'), b=get(b,'Signal');
end

% one of (a,b) must be an estruct else fallback to the usual linspace
if ~isa(a,'estruct') && ~isa(b,'estruct')
  v=linspace(a,b,n);
  return
end

% create constant objects from scalar input, using the other object as template
if ~isa(a, 'estruct') & isnumeric(a) & ~isempty(b)
  s=mean(a(:))*ones(size(get(b,'Signal'))); a=copyobj(b); set(a,'Signal', s);
elseif ~isa(b, 'estruct') & isnumeric(b) & ~isempty(a)
  s=mean(b(:))*ones(size(get(a,'Signal'))); b=copyobj(a); set(b,'Signal', s);
else
  [a,b] = intersect(a,b); % get intersection for operation (not needed when using copyobj)
end

xa = linspace(1,0,n);

v = zeros(estruct, n, 1);
for index=1:n
  this = a.*xa(index) + b.*(1-xa(index));
  this = setalias(this,'ratio',[ xa(index) 1-xa(index) ]);
  v(index) = this;
end


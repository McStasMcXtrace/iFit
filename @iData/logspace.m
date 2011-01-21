function v = logspace(a,b,n)
% v = logspace(a,b,n) : creates a logarithmically spaced vector of objects
%
%   @iData/logspace function to create a logarithmically spaced vector of iData objects
%   This corresponds to a 'morphing' from a to b in n steps.
%
% input:  a: object  (iData)
%         b: object  (iData)
%         n: number of steps, i.e. length of vector. Default is n=10 (integer)
% output: v: vector (iData array)
% ex:     b=logspace(a,b);
%
% Version: $Revision: 1.2 $
% See also iData, iData/max, iData/min, iData/colon, iData/linspace

if ~isa(a, 'iData') | ~isa(b,'iData')
  iData_private_error(mfilename, 'Operation requires 2 single iData objects as argument');
end

if nargin == 2
  n = [];
end
if isempty(n) | n <=0
  n=10;
end

if     isempty(a), a=0; 
elseif isempty(b), b=0; 
elseif isscalar(a) & isa(a, 'iData'), a=get(a,'Signal');
elseif isscalar(b) & isa(b, 'iData'), b=get(b,'Signal');
end

if ~isa(a, 'iData') & isscalar(a) & ~isempty(b)
  s=a*ones(size(get(b,'Signal'))); a=copyobj(b); set(a,'Signal', s);
elseif ~isa(b, 'iData') & isscalar(b) & ~isempty(a)
  s=b*ones(size(get(a,'Signal'))); b=copyobj(a); set(b,'Signal', s);
end


[a,b] = intersect(a,b);

xa = logspace(1,0,n);
xa = (xa-1); xa=xa/max(xa);

v = [];
for index=1:n
  c = a.*xa(index) + b.*(1-xa(index));
  v = [ v c ];
end


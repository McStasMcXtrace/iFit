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
% Version: $Revision: 1.1 $
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

[a,b] = intersect(a,b);

xa = logspace(1,0,n);

v = [];
for index=1:n
  c = a.*(xa(index)/10) + b.*((10-xa(index))/10);
  v = [ v c ];
end


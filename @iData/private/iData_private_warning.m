function iData_private_warning(a,b)
% for compatibility with Matlab < 6.5
if nargin == 1
  b=a;
  a='iData';
end
b = [ 'iData/' a ': ' b ];
a = [ 'iData:' a ];
try
  warning(a,sprintf(b));
catch
  warning(sprintf(b));
end

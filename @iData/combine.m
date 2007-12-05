function c = combine(a,b)
% c = combine(a,b) : combines iData objects
%
%   @iData/combine function to combine data sets
%     A fast notation for combine(a,b) is a\b
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array (iData)
% ex:     c=combine(a,b); or combine([ a b ])
%
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide
if nargin >= 1
  if ~isa(a, 'iData')
    iData_private_error(mfilename,['1st argument must be an iData object. Currently ' class(a) ]);
  elseif length(a) > 1, a=a(:); end
end
if nargin == 2
  if ~isa(a, 'iData')
    iData_private_error(mfilename,['2nd argument must be an iData object. Currently ' class(b) ]);
  elseif length(b) > 1, b=b(:); end
end
if nargin == 2
  c = combine([ a ; b ]); 
  return 
end

% now we should only handle a single argument
c = a(1);
if length(a) <= 1, return; end
for index=2:length(a)
  c = iData_private_binary(c, a(index), 'combine');
end



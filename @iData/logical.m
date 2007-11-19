function a=logical(a)
% d=logical(s) : convert iData into logical array
%
%   @iData/logical function to convert iData object into single floats
%
% input:  s: object or array (iData)
% output: v: value of the iData Signal (logical)
% ex:     'logical(iData(rand(10)))'
%  
% See also  iData/cell, iData/double, iData/struct, 
%           iData/char, iData/size

% EF 11/07/00 creation
% EF 23/09/07 iData implementation

if length(a) > 1
  b = {};
  for index=1:length(a(:))
    a{index} = iData_private_unary(a(index), op);
  end
  a = reshape(b, size(a));
  return
end

a = get(a, 'Signal');
a = logical(a);

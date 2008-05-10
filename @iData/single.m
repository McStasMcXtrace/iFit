function a=single(a)
% d = single(s) : convert iData into single floats
%
%   @iData/double function to convert iData object into single floats
%
% input:  s: object or array (iData)
% output: v: value of the iData Signal (single)
% ex:     'single(iData(rand(10)))'
%  
% See also  iData/cell, iData/double, iData/struct, 
%           iData/char, iData/size

% EF 11/07/00 creation
% EF 23/09/07 iData implementation

if length(a) > 1
  b = {};
  for index=1:length(a(:))
    b{index} = iData_private_unary(a(index), op);
  end
  a = reshape(b, size(a));
  return
end

a = get(a, 'Signal');
a = single(a);


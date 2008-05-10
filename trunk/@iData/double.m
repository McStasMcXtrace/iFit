function a=double(a)
% d = double(s) : convert iData into doubles
%
%   @iData/double function to convert iData object into doubles
%
% input:  s: object or array (iData)
% output: v: value of the iData Signal (double)
% ex:     'double(iData(rand(10)))'
%  
% See also  iData/cell, iData/single, iData/struct, 
%           iData/char, iData/size

% EF 11/07/00 creation
% EF 23/09/07 iData implementation

if length(a) > 1
  b = {};
  for index=1:length(a(:))
    b{index} = double(a(index));
  end
  b = reshape(b, size(a));
  return
end

a = get(a, 'Signal');
a = double(a);

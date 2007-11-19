function a=find(a)
% d=find(s) : find iData signal non zeros values
%
%   @iData/find function to find iData non zeros values
%
% input:  s: object or array (iData)
% output: v: value of the iData Signal (logical)
% ex:     'find(iData(rand(10)))'
%  
% See also  iData/find, iData

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
a = find(a);

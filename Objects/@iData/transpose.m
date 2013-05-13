function a = transpose(a)
% b = transpose(s) : non-conjugate transpose of iData object
%
%   @iData/transpose function to return the non-conjugate transpose of data sets
%   which corresponds to syntax: b = s.'
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=transpose(a);
%
% Version: $Revision: 1008 $
% See also iData, iData/transpose, iData/ctranspose, iData/setaxis, iData/getaxis

if numel(a) > 1
  a = builtin('transpose', a);
elseif ndims(a) <=2
  a = iData_private_unary(a, 'transpose');
else
  a = permute(a, [2 1]);
end


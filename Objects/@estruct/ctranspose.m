function a = ctranspose(a)
% b = ctranspose(s) : complex conjugate transpose of estruct object
%
%   @estruct/ctranspose function to return the complex conjugate transpose of data sets
%     which corresponds to syntax: b = s'
%   When the argument is an estruct array, the whole array is transposed. To
%     transpose each array element, use transpose(s) or b=s.'
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=ctranspose(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/transpose, estruct/ctranspose, estruct/setaxis, estruct/getaxis

if numel(a) > 1
  a = builtin('ctranspose', a);
elseif ndims(a) <=2
  a = unary(a, 'ctranspose');
else
  a = permute(a, [2 1]);
end

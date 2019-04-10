function a = transpose(a)
% b = transpose(s) : non-conjugate transpose of estruct object
%
%   @estruct/transpose function to return the non-conjugate transpose of data sets
%     which corresponds to syntax: b = s.'
%   When the argument is an estruct array, each element is transposed. To
%     transpose the array itself, use ctranspose(s) or b=s'
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=transpose(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/transpose, estruct/ctranspose, estruct/setaxis, estruct/getaxis

if numel(a) > 1
  a = unary(a, 'transpose'); % a = builtin('transpose', a);
elseif ndims(a) <=2
  a = unary(a, 'transpose');
else
  a = permute(a, [2 1]);
end


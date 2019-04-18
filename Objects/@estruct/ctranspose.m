function a = ctranspose(a)
%  '   Complex conjugate transpose.
%    X' is the complex conjugate transpose of X.
%    When the argument is an object array, the whole array is transposed. To
%    transpose each array element, use syntax: TRANSPOSE(s) or b=s.'
%
%    B = CTRANSPOSE(A) is called for the syntax A' (complex conjugate
%    transpose).
%
% Example: s=estruct(rand(3,6)); all(all(s' == s{0}'))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/transpose, estruct/ctranspose, estruct/setaxis, estruct/getaxis

if numel(a) > 1
  a = builtin('ctranspose', a);
elseif ndims(a) <=2
  a = unary(a, 'ctranspose');
else
  a = permute(a, [2 1]);
end
function a = transpose(a)
% .' Transpose.
%    X.' is the non-conjugate transpose.
%    When the argument is an array, each element is transposed. To
%    transpose the array itself, use CTRANSPOSE(s) or b=s'
%
%    B = TRANSPOSE(A) is called for the syntax A.'
%
% Example: s=estruct(rand(3,6)); all(size(transpose(s)) == [6 3])
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/transpose, estruct/ctranspose, estruct/setaxis, estruct/getaxis

if numel(a) > 1
  a = unary(a, 'transpose'); % a = builtin('transpose', a);
elseif ndims(a) <=2
  a = unary(a, 'transpose');
else
  a = permute(a, [2 1]);
end

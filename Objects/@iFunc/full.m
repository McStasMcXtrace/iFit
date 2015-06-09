function a = full(a)
% b = full(s) : Convert iFunc object storage to full matrix
%
%   @iFunc/full function to use full matrix storage.
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=full(a);
%
% Version: $Date$
% See also iFunc, iFunc/full, iFunc/sparse

a = iFunc_private_unary(a, 'full');

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end


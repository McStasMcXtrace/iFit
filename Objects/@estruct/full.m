function a = full(a)
% b = full(s) : Convert estruct object storage to full matrix
%
%   @estruct/full function to use full matrix storage, which stores
%   all elements in Signal, Error and Monitor (as opposed to sparse). 
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=full(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/pack, estruct/sparse

a = unary(a, 'full');

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end


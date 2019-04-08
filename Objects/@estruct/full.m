function a = full(a)
% b = full(s) : Convert iData object storage to full matrix
%
%   @iData/full function to use full matrix storage, which stores
%   all elements in Signal, Error and Monitor (as opposed to sparse). 
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=full(a);
%
% Version: $Date$ $Version$ $Author$
% See also iData, iData/pack, iData/sparse

a = iData_private_unary(a, 'full');

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end


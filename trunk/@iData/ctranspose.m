function a = ctranspose(a)
% b = ctranspose(s) : complex conjugate transpose of iData object
%
%   @iData/ctranspose function to return the complex conjugate transpose of data sets
%   which corresponds to syntax: b = s'
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=ctranspose(a);
%
% See also iData, iData/transpose, iData/ctranspose, iData/setaxis, iData/getaxis

a = iData_private_unary(a, 'ctranspose');


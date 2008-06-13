function c = power(a,b)
% c = power(a,b) : computes the power of iData objects
%
%   @iData/power function to compute the power of data sets
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array (iData)
% ex:     c=a.^2;
%
% Version: $Revision: 1.2 $
% See also iData, iData/times, iData/rdivide, iData/power

c = iData_private_binary(a, b, 'power');


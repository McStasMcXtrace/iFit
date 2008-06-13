function c = minus(a,b)
% c = minus(a,b) : computes the difference of iData objects
%
%   @iData/minus function to compute the difference of data sets
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array (iData)
% ex:     c=a-1;
%
% Version: $Revision: 1.2 $
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide

c = iData_private_binary(a, b, 'minus');


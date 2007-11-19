function c = plus(a,b)
% c = plus(a,b) : computes the sum of iData objects
%
%   @iData/plus function to compute the sum of data sets 
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array (iData)
% ex:     c=a+1;
%
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide

c = iData_private_binary(a, b, 'plus');


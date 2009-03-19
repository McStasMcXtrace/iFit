function c = ldivide(a,b)
% c = ldivide(a,b) : computes the left division ratio of iData objects
%
%   @iData/ldivide (.\) function to compute the left division ratio of data sets
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array (iData)
% ex:     c=a.\2; c=a.\b;
%
% Version: $Revision: 1.1 $
% See also iData, iData/minus, iData/plus, iData/times, iData/ldivide
if nargin ==1
	b=[];
end
c = iData_private_binary(a, b, 'ldivide');


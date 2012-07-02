function c = plus(a,b)
% c = plus(a,b) : computes the sum of iFunc objects
%
%   @iFunc/plus (+) function to compute the sum of functions
%
% input:  a: object or array (iFunc or numeric)
%         b: object or array (iFunc or numeric)
% output: c: object or array (iFunc)
% ex:     c=a+1;
%
% Version: $Revision: 1.1 $
% See also iFunc, iFunc/minus, iFunc/plus, iFunc/times, iFunc/rdivide

if nargin ==1
	b=[];
end
c = iFunc_private_binary(a, b, 'plus');


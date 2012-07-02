function c = rdivide(a,b)
% c = rdivide(a,b) : computes the ratio of iFunc objects
%
%   @iFunc/rdivide (./) function to compute the product of functions
%
% input:  a: object or array (iFunc or numeric)
%         b: object or array (iFunc or numeric)
% output: c: object or array (iFunc)
% ex:     c=lorz./gauss;
%
% Version: $Revision: 1.1 $
% See also iFunc, iFunc/minus, iFunc/plus, iFunc/times, iFunc/rdivide

if nargin ==1
	b=[];
end
c = iFunc_private_binary(a, b, 'rdivide');


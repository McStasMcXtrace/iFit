function c = minus(a,b)
% c = minus(a,b) : computes the difference of estruct objects
%
%   @estruct/minus (-) function to compute the difference of data sets
%
% input:  a: object or array (estruct or numeric)
%         b: object or array (estruct or numeric)
% output: c: object or array (estruct)
% ex:     c=a-1;
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide

if nargin ==1
	b=[];
end
c = binary(a, b, 'minus');


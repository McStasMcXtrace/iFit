function c = mrdivide(a,b)
% c = mrdivide(a,b) : computes the division of estruct objects
%
%   @estruct/mrdivide (/) function to compute the division of data sets=a/b
%
% input:  a: object or array (estruct or numeric)
%         b: object or array (estruct or numeric)
% output: c: object or array (estruct)
% ex:     c=a/2;
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide, estruct/power
if nargin ==1
	b=[];
end
c = binary(a, b, 'rdivide');



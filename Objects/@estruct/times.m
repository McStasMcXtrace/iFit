function c = times(a,b)
% c = times(a,b) : computes the product of estruct objects
%
%   @estruct/times (*) function to compute the product of data sets=a.*b
%     the square of a single estruct object should rather be computed 
%     using the power law.
%
% input:  a: object or array (estruct or numeric)
%         b: object or array (estruct or numeric)
% output: c: object or array (estruct)
% ex:     c=a.*2;
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide, estruct/power
if nargin ==1
	b=[];
end

c = binary(a, b, 'times');


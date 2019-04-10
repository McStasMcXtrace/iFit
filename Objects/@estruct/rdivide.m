function c = rdivide(a,b)
% c = rdivide(a,b) : computes the ratio of estruct objects
%
%   @estruct/rdivide (./) function to compute the ratio of data sets
%
% input:  a: object or array (estruct or numeric)
%         b: object or array (estruct or numeric)
% output: c: object or array (estruct)
% ex:     c=a./2; c=a./b;
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide
if nargin ==1
	b=[];
end
c = binary(a, b, 'rdivide');


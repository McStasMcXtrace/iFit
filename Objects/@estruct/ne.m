function c = ne(a,b)
% c = ne(a,b) : not-equal comparison of estruct objects
%
%   @estruct/ne (~=) comparison operator
%     when comparing two estruct objects, the monitor weighting is applied.
%
% input:  a: object or array (estruct or numeric)
%         b: object or array (estruct or numeric)
% output: c: object or array which Signal is the comparison result (estruct)
% ex:     c= (a~=1); d=find(a~=b);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/find, estruct/gt, estruct/lt, estruct/ge, estruct/le, estruct/ne, estruct/eq

if nargin ==1
	b=[];
end
c = binary(a, b, 'ne');


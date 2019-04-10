function c = plus(a,b)
% c = plus(a,b) : computes the sum of estruct objects
%
%   @estruct/plus (+) function to compute the sum of data sets 
%
%     the sum is defined as the normalised sum:
%       (M1+M2)*(S1/M1+S2/M2) over monitor(M1+M2)
%     where S1,M1 and S2,M2 and the Signal and Monitor of the two objects.
%
%     To get the 'conventional' sum which is the sum S1+S2, set one of the
%     monitors to 0, e.g.:
%       a.Monitor=0;
%
% input:  a: object or array (estruct or numeric)
%         b: object or array (estruct or numeric)
% output: c: object or array (estruct)
% ex:     c=a+1;
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide, estruct/combine

if nargin ==1
	b=[];
end
c = binary(a, b, 'plus');


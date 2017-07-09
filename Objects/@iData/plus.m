function c = plus(a,b)
% c = plus(a,b) : computes the sum of iData objects
%
%   @iData/plus (+) function to compute the sum of data sets 
%
%     the sum is defined as the normalised sum:
%       (M1+M2)*(S1/M1+S2/M2) over monitor(M1+M2)
%     where S1,M1 and S2,M2 and the Signal and Monitor of the two objects.
%
%     To get the 'conventional' sum which is the sum S1+S2, set one of the
%     monitors to 0, e.g.:
%       a.Monitor=0;
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array (iData)
% ex:     c=a+1;
%
% Version: $Date$
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide, iData/combine

if nargin ==1
	b=[];
end
c = iData_private_binary(a, b, 'plus');


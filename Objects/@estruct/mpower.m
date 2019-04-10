function c = mpower(a,b)
% c = mpower(a,b) : computes the matrix-power of estruct objects
%
%   @estruct/mpower function to compute the power of data sets
%
% input:  a: object or array (estruct or numeric)
%         b: object or array (estruct or numeric)
% output: c: object or array (estruct)
% ex:     c=a^2;
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/times, estruct/rdivide, estruct/power

if nargin == 1,
  b = a;
end
c = binary(a, b, 'mpower');


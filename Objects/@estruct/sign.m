function a = sign(a)
% b = sign(s) : sign of estruct object
%
%   @estruct/sign function to return the sign of data sets
%   This function computes the sign of the object 's', i.e
%   -1 for negative values, 0 for null, and +1 for positive values.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=sign(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse

a = unary(a, 'sign');


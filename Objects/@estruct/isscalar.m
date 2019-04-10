function a = isscalar(a)
% b = isscalar(s) : True for scalar estruct objects
%
%   @estruct/isscalar function to return true if data set is a single scalar element
%   of 's', i.e. that size is 1x1
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=isscalar(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse

a = unary(a, 'isscalar');


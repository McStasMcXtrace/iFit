function a = issparse(a)
% b = issparse(s) : True for sparse matrix estruct object
%
%   @estruct/issparse function to return true if data set is a vector
%   of 's', i.e. that size is 1xN or Nx1
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=issparse(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse

a = unary(a, 'issparse');
a = uint8(a);

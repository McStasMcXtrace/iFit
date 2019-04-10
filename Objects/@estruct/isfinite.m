function a = isfinite(a)
% b = isfinite(s) : True for finite estruct object elements
%
%   @estruct/isfinite function to return true for finite elements
%   of 's', i.e. that are not NaN, Inf or -Inf.
%
%   To remove nan's and inf's values use: fill(s)
%
% input:  s: object or array (estruct)
% output: b: array (int)
% ex:     b=isfinite(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse, estruct/fill

a = unary(a, 'isfinite');


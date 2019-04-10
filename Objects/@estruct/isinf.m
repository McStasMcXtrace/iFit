function a = isinf(a)
% b = isinf(s) : True for infinite estruct object elements
%
%   @estruct/isinf function to return true for infinite elements
%   of 's', i.e. that are +Inf or -Inf 
%
%   To remove nan's and inf's values use: fill(s)
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=isinf(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse, estruct/fill
a = unary(a, 'isinf');


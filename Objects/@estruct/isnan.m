function a = isnan(a)
% b = isnan(s) : True for NaN estruct object elements
%
%   @estruct/isnan function to return true for NaN elements
%   of 's', i.e. that are NaN ('not a number')
%
%   To remove nan's and inf's values use: fill(s)
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=isnan(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/sign, estruct/isreal, estruct/isfinite, estruct/isnan,
%          estruct/isinf, estruct/isfloat, estruct/isinterger,
%          estruct/isnumeric, estruct/islogical, estruct/isscalar, 
%          estruct/isvector, estruct/issparse, estruct/fill

a = unary(a, 'isnan');

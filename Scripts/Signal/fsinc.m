function s=fsinc(x)
% s=sinc(x)
%   returns an array, y, whose elements are the sinc of the elements 
%   of the input, x. y=sin(x)/x is the same size as x.
%
% (c) E.Farhi, ILL. License: EUPL.
if ndims(x) > 1, return; end

index    = find(x ~= 0 & isfinite(x));
s        = zeros(size(x));
s(index) = sin(x(index)) ./ x(index);

% special case for x=0
s(x == 0) = 1;

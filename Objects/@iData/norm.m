function v = norm(a, varargin)
% NORM   Object norm.
%   Computes the norm of the object Signal. The default is to use the 
%   2-norm, defined as sqrt(sum( |a|^2 ))
%
%     NORM(V,P)    = sum(abs(V).^P)^(1/P).
%     NORM(V)      = norm(V,2).
%     NORM(V,inf)  = max(abs(V)).
%     NORM(V,-inf) = min(abs(V)).
%
% Example: X=iData([0 1 2 3]); round(norm(X))==4
% Version: $Date$ $Version$ $Author$
% See also iData, iData/sum, iData/trapz, norm

v = unary(a, 'norm', varargin{:});


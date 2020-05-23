function [s,sigma] = prod(a,dim)
% PROD Product of elements.
%    S = PROD(X) is the product of the Signal of the object X. If
%    X is a N-D object, PROD(X) operates along the first
%    dimension.
%    If X is floating point, that is double or single, S is
%    multiplied natively, that is in the same class as X,
%    and S has the same class as X. If X is not floating point,
%    S is multiplied in double and S has class double.
% 
%    PROD(X,DIM) works along the dimension DIM. 
%
%    [S,sigma] = PROD(X, 0) does the same as above, but returns the total 
%    product per object. In this case, a second output argument holds the error
%    bar.
%
% Example: s=estruct(-10:10); prod(s,0) == 0
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plus, estruct/prod, estruct/cumprod, estruct/mean

if nargin < 2, dim=1; end

[s,sigma] = private_sumtrapzproj(a,dim, 'prod');

if nargin > 1 && isequal(dim,0)
  if iscell(s)
    s = cell2mat(s);
  end
end

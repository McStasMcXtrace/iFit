function a = full(a)
%  FULL   Convert sparse object to full matrix.
%    A = FULL(X) converts a sparse object S to full storage
%    organization.  If X is a full object, it is left unchanged.
%    FULL affects Signal, Error and Monitor in the object.
%
%    If A is full, issparse(A) returns 0.
%
% Example: S = iData(sparse(rand(200,200) < 2/3)); F=full(S); ~issparse(F) & issparse(S)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/pack, iData/sparse, iData/pack

a = unary(a, 'full');
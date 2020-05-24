function a = sparse(a)
%  SPARSE Create sparse matrix.
%    S = SPARSE(X) converts a sparse or full object to sparse form by
%    squeezing out any zero elements. This may be usefull for event
%    based storage where most events are zeros. Use iData/full to
%    revert to full matrix storage. The Signal, Error and Monitor
%    are affected.
%
% Example: S = iData(sparse(rand(200,200) < 2/3)); issparse(S)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/full, iData/pack, iData/issparse

a = unary(a, 'sparse');
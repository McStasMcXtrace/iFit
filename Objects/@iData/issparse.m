function v = issparse(a)
%  ISSPARSE True for sparse matrix.
%    ISSPARSE(S) is 1 if the storage class of S Signal is sparse
%    and 0 otherwise.
%
% Example: s=iData(sparse(eye(5))); issparse(s)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/full, iData/sparse

v = unary(a, 'issparse');
if iscell(v), v=cell2mat(v); end
v = logical(v);
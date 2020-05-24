function v = islogical(a)
%  ISLOGICAL True for logical object.
%    ISLOGICAL(X) returns true if X is a logical object Signal and false otherwise.
%    Logical arrays must be used to perform logical 0-1 indexing.
%
% Example: s=iData([true false]); islogical(s)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/isnumeric, logical

v = unary(a, 'islogical');
if iscell(v), v = cell2mat(v); end
v = logical(v);
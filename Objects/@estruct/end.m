function b = end(s,k,n)
% END Last index in an indexing expression.
%   END(A,K,N) is called for indexing expressions involving the object A
%   when END is part of the K-th index out of N indices.
%
% Version: $Date$ $Version$ $Author$
% See also: estruct

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

if numel(s) > 1
  if n == 1, b=numel(s); else b=size(s,k); end
  return
end

if n > length(size(s))
  error([ mfilename ': input object ' inputname(1) ' ' b.Tag ' ' b.name ...
    ' has a size [' num2str(size(s)) '] but the dimension ' n ' is requested.' ]);
end

S = subsref(s,struct('type','.','subs','Signal'));

if n == 1 && ndims(S) > 1
  % special case object(end) always return the real last element
  b = prod(size(S));
else
  % last element along dimension
  b = size(S,k);
end
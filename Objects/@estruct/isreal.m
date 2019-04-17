function v = isreal(a)
%  ISREAL True for real object.
%    ISREAL(X) returns 1 if X does not have an imaginary part
%    and 0 otherwise.
%
%    ~ISREAL(X) detects objects that have an imaginary part even if
%    it is all zero.
%    ~ANY(IMAG(X),0) detects strictly real objects, whether X has
%    an all zero imaginary part allocated or not.
%
% Example: s=estruct(rand(5)+i*rand(5)); ~isreal(s)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/real, estruct/imag, complex

v = unary(a, 'isreal');
if iscell(v), v = cell2mat(v); end
v = logical(v);
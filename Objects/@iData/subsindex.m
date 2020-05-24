function y = subsindex(s)
% SUBSINDEX Subscript index.
%   I = SUBSINDEX(A) is called for the syntax 'X(A)' when A is an
%   object.  SUBSINDEX must return the value of the object as a
%   zero-based integer index (I must contain integer values in the
%   range 0 to prod(size(X))-1).  SUBSINDEX is called by the default
%   SUBSREF and SUBSASGN functions and you may call it yourself if you
%   overload these functions.
%
%    SUBSINDEX is invoked separately on all the subscripts in an
%    expression such as X(A,B).
%
% Example: X=iData(5:10); a=iData(-5:2:50); b=a(X); all(b == 1:2:11)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/subsasgn, iData/subsref

if numel(s) > 1
   s = s(1);
end


if ~islogical(s)
  y = uint32(s);
  if any(y(:) <= 0)
    error([ mfilename ': Invalid index values from object ' s.Tag ' ' s.Name '. Some indices are <= 0.' ]);
  end
else
  y = logical(s);
end
y = y(:)-1; % between 0 and max-1.

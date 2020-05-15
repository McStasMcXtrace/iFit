function b = not(a)
% ~   Logical NOT.
%   ~A performs a logical NOT of input object A, and returns an object
%    containing elements set to either logical 1 (TRUE) or logical 0 (FALSE).
%    An element of the output object is set to 1 if A Signal contains a zero value
%    element at that same array location.  Otherwise, that element is set to
%    0.
%
%    B = NOT(A) is called for the syntax '~A'.
%
% Example: X=estruct([1 0 1 0]); all(not(X)==[0 1 0 1])
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/or, estruct/and

b = unary(a, 'not');

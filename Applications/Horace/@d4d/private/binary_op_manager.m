% Implement binary arithmetic operations for objects containing a double array.
%
%   if w1, w2 are objects of the same size:
%       - the operation is performed element-by-element
%
%   if one of w1 or w2 is double:
%       - if a scalar, apply to each element of the object double array
%       - if an array of the same size as the object double array, apply
%        element by element
%
%   w1, w2 can be arrays:
%       - if objects have same array sizes, then add element-by-element
%       - if an (n+m)-dimensional array, the inner n dimensions will be
%        combined element by element with the object double array (where
%        n is the dimensionality of the object double array), and the
%        outer m dimensions must match the array size of the array of objects
%
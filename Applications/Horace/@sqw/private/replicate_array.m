% Replicate array elements according to list of repeat indicies
%
%   >> vout = replicate_array (v, n)
%
%   v       Array of values
%   n       List of number of times to replicate each value
%
%   vout    Output array: column vector
%               vout=[v(1)*ones(1:n(1)), v(2)*ones(1:n(2), ...)]'
%
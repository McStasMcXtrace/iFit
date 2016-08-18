% Implement w1 - w2 for objects
%
%   >> w = w1 - w2
%
%   if w1, w2 are objects of the same size:
%       - the operation is performed element-by-element
%
%   if one of w1 or w2 is numeric:
%       - if a scalar, apply to each element of the object numeric array
%       - if an array of the same size as the object numeric array, apply
%        element by element
%
%   w1, w2 can be arrays:
%       - if objects have same array sizes, then add element-by-element
%       - if an (n+m)-dimensional array, the inner n dimensions will be
%        combined element by element with the object numeric array (where
%        n is the dimensionality of the object numeric array), and the
%        outer m dimensions must match the array size of the array of objects
%%   Overloaded methods:
%      sqw/minus
%      testsigvar/minus
%      sigvar/minus
%      IX_dataset_3d/minus
%      IX_dataset_2d/minus
%      IX_dataset_1d/minus
%      sqw/minus
%      sigvar/minus
%      d4d/minus
%      d3d/minus
%      d2d/minus
%      d1d/minus
%      d0d/minus
%      codistributed/minus
%      gpuArray/minus
%      timeseries/minus
%%   Reference page in Help browser
%      doc minus
%
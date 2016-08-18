% Implement w1 \ w2 for objects
%
%   >> w = w1 \ w2
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
%      sqw/mldivide
%      testsigvar/mldivide
%      sigvar/mldivide
%      IX_dataset_3d/mldivide
%      IX_dataset_2d/mldivide
%      IX_dataset_1d/mldivide
%      sqw/mldivide
%      sigvar/mldivide
%      d4d/mldivide
%      d3d/mldivide
%      d2d/mldivide
%      d1d/mldivide
%      d0d/mldivide
%      codistributed/mldivide
%      gpuArray/mldivide
%      timeseries/mldivide
%%   Reference page in Help browser
%      doc mldivide
%
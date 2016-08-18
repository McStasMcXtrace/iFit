% Get limits of the data in an n-dimensional dataset, that is, find the
% coordinates along each of the axes of the smallest cuboid that contains
% bins with non-zero values of contributing pixels.
%
% Syntax:
%   >> [val, n] = data_bin_limits (din)
%
% Input:
% ------
%   din     Input dataset structure
%
% Output:
% -------
%   val     (2 x ndim) array, where ndim = dimension of dataset,containing
%           the lower and upper limits of the bin boundaries of the dataset.
%           isempty(val)=1 if there is no data in the dataset
%   
%   n       (2 x ndim) array containing the lower and upper indices of the 
%           elements along each axis
%           isempty(n)=true if there is no data in the dataset
%
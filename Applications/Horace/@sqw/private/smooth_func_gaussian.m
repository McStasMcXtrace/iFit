% Fills a matrix with a Gaussian function in multiple dimensions
% Output matrix normalised so sum of elements = 1
%
% Syntax:
%   >> c = smooth_func_gaussian(width)
%
%   width   Vector containing full-width half-height (FWHH) in pixels of Gaussian 
%           function along each dimension.
%           The convolution matrix extends to 2% of the peak intensity 
%
%   c       Convolution array
%
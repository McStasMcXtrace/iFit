% Permute the pixel information array elements according to list of permutations of the plot axes
% Arrays pix and npix are assumed to be compatible, and order too - not all checks are done !
%
%   >> pixout = permute_pix_array (pix, npix, order)
%
%   pix     Array of pixel information values (9 x npixtot)
%   npix    List of number of pixels contributing to each bin.
%   order   Rearrangement of axes, foolowing usual convention of matlab intrinsic PERMUTE
%          e.g. if npix is 3D with size [13,17,22], then 
%          order [2,3,1] rearranges pix so that it corresponds to a new npix with size [17,22,13]
%
%   pixout  Output array (9 x npixtot)
%
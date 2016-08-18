% Rewmove blocks from an array
%
%   >> vout = compress_array (v, npix, mask)
%
%   v       Two dimensional array of values, v(:,ntot) where
%              ntot=sum(n)
%   n       List of number of elements along outer dimension of v
%          corresponding to one block
%   mask    Block indicies to remove. Must have size(mask)==size(n)
%
%   vout    Output array, v(:,ntot_compress), where
%               ntot_compress = sum(n(~mask))
%
function s = read_image(filename)
% read_image Wrapper to imfinfo/imread which reconstructs the image structure
%   s = read_image(filename)
%
% This function can read (imformats): 
%   fits gif hdf jpeg pbm png tiff ...
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_fits, imformats

s=[];
if nargin == 0, return; end

s       = imfinfo(filename);
s.image = imread(filename);
if exist('exifread') == 2
    try
    s.EXIF = exifread(filename);
    end
end


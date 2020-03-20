function s = read_image(filename)
% read_image Wrapper to imfinfo/imread which reconstructs the image structure
%   s = read_image(filename)
%
% This function can read (imformats): 
%   fits gif hdf jpeg pbm png tiff ...
%
% Input:  filename: image file (string)
% output: structure
% Example: y=read_image(fullfile(ifitpath, 'Data','Ag_3_a.hdf4')); isstruct(y)
%
% 
% See also: read_fits, imformats

s=[];


if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    % we get the imformats extensions
    imf = imformats;
    imfext = {};
    for index=1:numel(imf)
      imfext = [ imfext imf(index).ext{:} ];
    end
    
    s.name           = 'Image/Picture';
    s.extension      = imfext;
    s.method         = mfilename;
    return
end

try
  s       = imfinfo(filename);
catch
  s = [];
  return;
end
s.image = imread(filename);
if exist('exifread') == 2
    warning('off','MATLAB:exifread:DeprecatedFunction');
    try
    s.EXIF = exifread(filename);
    end
end


function s = read_fits(filename)
% read_fits Wrapper to fitsinfo/fitsread which reconstructs the FITS structure
%   s = read_fits(filename)
%
% Input:  filename: FITS image file (string)
% output: structure
% Example: y=read_fits(fullfile(ifitpath, 'Data','30dor.fits')); isstruct(y)
%
% 
% See also: read_edf, read_adsc, read_edf, read_sif, read_mar, read_spe, read_cbf, read_image, read_hbin

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    s.name     ='FITS';
    s.method   =mfilename;
    s.extension={'fits','fts'};
    s.patterns ={'SIMPLE','BITPIX','NAXIS','END'};
    return
end

try
  s      = fitsinfo(filename);
catch
  s=[]; return;
end
% we read all FITS 'extnames'
s.data = [];
for extname={'primary','table','bintable','image','unknown'}
  try
    s.data.(extname{1}) = fitsread(filename, extname{1});
  end
end


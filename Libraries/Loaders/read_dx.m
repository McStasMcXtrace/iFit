function s = read_dx(filename)
% READ_DX Read an Electronic Chemistry Exchange format file, used e.g. in
% infrared, NMR, and mass spectroscopy.
%   s = read_dx(filename)
%
% http://www.jcamp-dx.org/
%
% See also: read_jeol, read_bruker

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    dx.name            ='Electronic Chemistry Exchange format (.dx)';
    dx.extension       ='dx';
    dx.method          =mfilename;
    % dx.postprocess     ={'opencdf','load_nmoldyn'};
    s = dx;
    return
end

s = tsread(filename);

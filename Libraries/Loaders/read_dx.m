function s = read_dx(filename)
% READ_DX Read an Electronic Chemistry Exchange format file, used e.g. in
% infrared, NMR, and mass spectroscopy.
% These are also referred as Joint Committee on Atomic and Molecular Physical data.
%   s = read_dx(filename)
%
% Example: y=read_dx(fullfile(ifitpath, 'Data','2Methyl1Propanol.jdx')); isstruct(y)
%
% References:
%   http://www.jcamp-dx.org/
%   http://www.chm.bris.ac.uk/~paulmay/temp/pcc/jcamp.htm
%
% See also: read_jeol, read_bruker, read_varian, read_opus

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    dx.name            ='Electronic Chemistry Exchange/JCAMP format (.dx)';
    dx.extension       = {'dx','jdx'};
    dx.method          = mfilename;
    % dx.postprocess     ={'opencdf','load_nmoldyn'};
    s = dx;
    return
end

s = tsread(filename);

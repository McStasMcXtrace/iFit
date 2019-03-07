function s = read_ply(filename, varargin)
% data=read_ply(filename, options, ...) Read PLY 3D stereo-lithography CAD data
%
% read_ply Read PLY 3D stereo-lithography CAD data
%
%   Reference: PLY format at https://en.wikipedia.org/wiki/PLY_%28file_format%29
%
% Input:  filename: PLY Data text file (string)
% output: structure
% Example: y=read_ply(fullfile(ifitpath, 'Data','socket.ply')); isstruct(y)
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_stl, read_off, read_anytext

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    PLY_ascii.name      ='PLY 3D ascii';
    PLY_ascii.method    =mfilename;
    PLY_ascii.options   ='--fast --binary --headers --comment=NULL --silent ';
    PLY_ascii.extension ='ply';
    PLY_ascii.patterns  ={'ply','format ascii','element','end_header'};
    PLY_ascii.postprocess='openply';
    
    s = PLY_ascii;
    return
end

% now call read_anytext with given options

if isempty(varargin)
  varargin = { '--fast --binary --headers --comment=NULL --silent ' };
end
s       = read_anytext(filename, varargin{:});

end


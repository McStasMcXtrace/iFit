function s = read_xye(filename, varargin)
% data=read_xye(filename, options, ...) Read Simple x/y/e column file
%
% read_xye Read Simple x/y/e column file
%
% Input:  filename: Simple x/y/e column file (string)
% output: structure
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_anytext

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    xye.name            ='Simple x/y/e column file';
    xye.extension       ={'xye', 'xy'};
    xye.method          =mfilename;
    xye.patterns        ={};
    xye.postprocess     = 'load_xyen';
    
    s = xye;
    return
end

% now call read_anytext with given options

s       = read_anytext(filename, varargin{:});

end


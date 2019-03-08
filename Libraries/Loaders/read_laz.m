function s = read_laz(filename, varargin)
% data=read_laz(filename, options, ...) Read McCode LAZ/LAU structure file
%
% read_laz read a McCode LAZ (powder Lazy/Pulvrx) or LAU (single crystal)
%   structure file.
%
% Input:  filename: McCode LAZ/LAU text file (string)
% output: structure
% Example: y=read_laz(fullfile(ifitpath, 'Data','Al.lau')); isstruct(y)
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_mccode, read_anytext

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    mcstas_powder.name       ='McCode powder table (LAZ/LAU)';
    mcstas_powder.patterns   ={'lattice_a','column_'};
    mcstas_powder.options    ='--fast --binary  --headers --comment=NULL --silent ';
    mcstas_powder.method     =mfilename;
    mcstas_powder.postprocess='openlaz';
    mcstas_powder.extension  ={'laz','lau'};
    
    s = mcstas_powder;
    return
end

% now call read_anytext with given options

if isempty(varargin)
  varargin = { '--fast --binary  --headers --comment=NULL --silent ' };
end
s       = read_anytext(filename, varargin{:});

end


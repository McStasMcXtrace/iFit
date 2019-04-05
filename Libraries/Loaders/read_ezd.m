function s = read_ezd(filename, varargin)
% data=read_ezd(filename, options, ...) Read EZD Electron density maps
%
% read_ezd Read EZD Electron density maps
%
%   Reference: EZD format at http://www.msg.ucsf.edu/local/programs/ono/manuals/ofaq/Q.803.html
%
% Input:  filename: EZD Data text file (string)
% output: structure
% Example: y=read_ezd(fullfile(ifitpath, 'Data','GroEl.ezd')); isstruct(y)
%
% $
% See also: read_mrc, read_anytext

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    EZD.name            ='EZD electronic density map';
    EZD.method          =mfilename;
    EZD.options         ='--fortran --headers --binary --fast --catenate --comment=NULL --silent';
    EZD.extension       ='ezd';
    EZD.patterns        ={'EZD_MAP','CELL','EXTENT'};
    EZD.postprocess     ='this.Data.MAP = reshape(this.Data.MAP, this.Data.EXTENT);';
    
    s = EZD;
    return
end

% now call read_anytext with given options

if isempty(varargin)
  varargin = { '--fortran --headers --binary --fast --catenate --comment=NULL --silent' };
end
s       = read_anytext(filename, varargin{:});

end


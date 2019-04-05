function s = read_off(filename, varargin)
% data=read_off(filename, options, ...) Read OFF 3D stereo-lithography CAD data
%
% read_off Read OFF 3D stereo-lithography CAD data
%
%   Reference: OFF format at http://www.geomview.org/docs/html/OFF.html
%
% Input:  filename: OFF Data text file (string)
% output: structure
% Example: y=read_off(fullfile(ifitpath, 'Data','socket.off')); isstruct(y)
%
% $
% See also: read_stl, read_ply, read_anytext

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    OFF_ascii.name      ='OFF 3D stereo-lithography CAD ascii';
    OFF_ascii.method    =mfilename;
    OFF_ascii.options   ='--fast --binary --headers --comment=NULL --metadata=OFF --silent ';
    OFF_ascii.extension ='off';
    OFF_ascii.patterns  ={'\<OFF\>'};
    OFF_ascii.postprocess='openoff';
    
    s = OFF_ascii;
    return
end

% now call read_anytext with given options

if isempty(varargin)
  varargin = { '--fast --binary --headers --comment=NULL --metadata=OFF --silent ' };
end
s       = read_anytext(filename, varargin{:});

end


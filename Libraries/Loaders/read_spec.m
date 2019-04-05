function s = read_spec(filename, varargin)
% data=read_spec(filename, options, ...) Read synchrotron/x-ray SPEC data
%
% read_spec Read synchrotron/x-ray SPEC data format
%
%   Reference: SPEC at Certified Scientific Software https://certif.com/
%
% Input:  filename: SPEC Data text file (string)
% output: structure
% Example: y=read_spec(fullfile(ifitpath, 'Data','SPEC.dat')); isstruct(y)
%
% 
% See also: read_llb_tas, read_anytext, read_ill

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    spec.name           ='SPEC';
    spec.patterns       ={'#F','#D','#S'};
    spec.options        ='--fast --binary --headers --metadata=''#S '' --comment=NULL --silent ';
    spec.method         =mfilename;
    spec.extension      ={'spc','spec','dat'};
    
    s = spec;
    return
end

% now call read_anytext with given options

if isempty(varargin)
  varargin = { '--fast --binary --headers --metadata=''#S '' --comment=NULL --silent ' };
end
s       = read_anytext(filename, varargin{:});

end


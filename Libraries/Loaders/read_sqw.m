function s = read_sqw(filename, varargin)
% data=read_sqw(filename, options, ...) Read SQW from McCode or ISIS Horace
%
% read_sqw read a McCode SQW (inelastic 2D/4D S(q,w)) or ISIS Horace ToF file
%
% Input:  filename: SQW McCode/Horace file (string)
% output: structure/object
% Example: y=read_sqw(fullfile(ifitpath, 'Data','Horace_sqw_1d.sqw')); isa(y, 'sqw')
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_mccode, read_anytext

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    ISIS_SQW.name       = 'ISIS Horace SQW';
    ISIS_SQW.method     = 'sqw';
    ISIS_SQW.extension  = 'sqw';
    ISIS_SQW.postprocess= 'this.Data = struct(this.Data); this = setaxis(this,0);';
    
    mcstas_sqw.name       ='McCode Sqw table';
    mcstas_sqw.patterns   ={'Sqw data file'};
    mcstas_sqw.options    ='--fast --binary  --headers --comment=NULL --silent --metadata=title ';
    mcstas_sqw.method     =mfilename;
    mcstas_sqw.postprocess='opensqw';
    mcstas_sqw.extension  ={'sqw','sqw4'};
    
    s = { mcstas_sqw, ISIS_SQW };
    return
end

% is the file binary ?
isbinary= 0;
% read start of file
[fid,message] = fopen(filename, 'r');
if fid == -1, return; end
file_start    = fread(fid, 1000, 'uint8=>char')';

% check if this is a text file
if length(find(file_start >= 32 & file_start < 127))/length(file_start) < 0.4
  isbinary = 1; % less than 90% of printable characters
end
fclose(fid);

% now call sqw for pure binary files (ISIS SQW), or read_anytext with given options
if isbinary && exist('sqw') && isa(sqw, 'sqw')
  try
    s = sqw(filename);
    return
  end
end

% else try a text reader
if isempty(varargin)
  varargin = { '--fast --binary  --headers --comment=NULL --silent --metadata=title ' };
end
s       = read_anytext(filename, varargin{:});




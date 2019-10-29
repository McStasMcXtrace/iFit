function s = read_pos(filename)
% read_pos Read Atom Probe Tomography defined written by e.g. Cameca LEAP instruments.
%   s = read_pos(filename)
%
% A POS file is a simple binary file (big endians) consisting of 4-byte 
%   IEEE float32s representing the x, y, and z position (in nm), and 
%   mass/charge (in amu) of every atom in a dataset. It is used by many 
%   atom probe groups and software packages and is the current de facto 
%   exchange method. The format is detailed in 
%   https://github.com/oscarbranson/apt-tools/blob/master/file-format-info.pdf
%
% Input:  filename: POS file (string)
% output: structure
% Example: y=read_pos(fullfile(ifitpath, 'Data','test.pos')); isstruct(y)
%
% 
% See also: read_nii, read_mrc, read_analyze

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    s.name            ='Atom Probe Tomography Exchange format (.pos)';
    s.extension       ='pos';
    s.method          =mfilename;
    return
end

% read the file
[fid,mess] = fopen(filename);
if fid == -1, error(mess); end

A=fread(fid, inf, 'float32','ieee-be');
fclose(fid);

% create the output data
s.x=A(1:4:end); s.y=A(2:4:end); s.z=A(3:4:end); s.mn=A(4:4:end);




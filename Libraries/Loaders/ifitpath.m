function d = ifitpath
% ifitpath iFit library location
%
% returns the iFit installation path (string)
%
% Example: d=ifitpath;
%
% Version: $Date$
% (c) E.Farhi, ILL. License: EUPL.

d = [ fileparts(which('iData/version')) filesep '..' filesep '..' filesep ];



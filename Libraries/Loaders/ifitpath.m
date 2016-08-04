function d = ifitpath
% ifitpath iFit library location
%
% Version: $Date$
% (c) E.Farhi, ILL. License: EUPL.

d = [ fileparts(which('iData/version')) filesep '..' filesep '..' filesep ];



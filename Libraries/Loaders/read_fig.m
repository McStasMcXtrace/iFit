function s = read_fig(filename)
% read_fig Wrapper to directly read Matlab Figures
%  s = read_fig(filename)
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_idl, read_lvm, read_tdms, read_igor

f       = openfig(filename, 'new','invisible');
s.Handle=f;


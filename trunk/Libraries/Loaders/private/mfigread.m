function s = mfigread(filename)
% mfigread Wrapper to directly read Matlab Figures

f       = openfig(filename, 'new','invisible');
s.Handle=f;


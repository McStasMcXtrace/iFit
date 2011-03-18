function s = mfigread(filename)
% mfigread Wrapper to directly read Matlab Figures

f      = openfig(filename, 'invisible');
s      = iData(f);
close(f);


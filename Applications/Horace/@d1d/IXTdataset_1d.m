function wout = IXTdataset_1d (w)
% Convert d1d object into IXTdataset_1d
%
%   >> wout = IXTdataset_1d (w)

% Original author: T.G.Perring
%
% $Revision: 301 $ ($Date: 2009-11-03 20:52:59 +0000 (Tue, 03 Nov 2009) $)

wout=IXTdataset_1d(sqw(w));

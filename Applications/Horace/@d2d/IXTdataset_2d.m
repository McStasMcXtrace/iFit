function wout = IXTdataset_2d (w)
% Convert d2d object into IXTdataset_2d
%
%   >> wout = IXTdataset_2d (w)

% Original author: T.G.Perring
%
% $Revision: 301 $ ($Date: 2009-11-03 20:52:59 +0000 (Tue, 03 Nov 2009) $)

wout=IXTdataset_2d(sqw(w));

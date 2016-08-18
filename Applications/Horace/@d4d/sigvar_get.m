function [s,var,mask_null] = sigvar_get (w)
% Get signal and variance from object, and a logical array of which values to ignore
% 
%   >> [s,var,mask_null] = sigvar_get (w)

% Original author: T.G.Perring
%
% $Revision: 301 $ ($Date: 2009-11-03 20:52:59 +0000 (Tue, 03 Nov 2009) $)

s = w.s;
var = w.e;
mask_null = logical(w.npix);

function b = version(a)
% v = version(iData): iData class version
%
%   @iData/version: version of the iData class library
%
% Version: $Revision: 1.10 $

% EF 23/09/07 iData impementation

vers = '@IFIT_VERSION@';
date = '@IFIT_DATE@';
auth = 'E.Farhi, P. Willendrup and Y.Debab, (c) ILL/DS/CS <farhi@ill.eu> GPL2';
b = [ 'Matlab/iData version ' vers ' (' date ') by ' auth ];



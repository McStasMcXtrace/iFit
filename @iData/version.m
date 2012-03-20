function b = version(a,long_request)
% v = version(iData): iData class version
%
%   @iData/version: version of the iData class library
%
% Version: $Revision: 1.14 $

% EF 23/09/07 iData impementation

vers = '@IFIT_VERSION@';
date = '@IFIT_DATE@';
auth = 'E.Farhi, P. Willendrup and Y.Debab, (c) ILL/DS/CS <farhi@ill.eu> GPL2';
contrib = 'Eric Ludlam, Felix Morsdorf, Joe Conti, Douglas M. Schwarz, Alexandros Leontitsis, F. Sigworth, Argimiro R. Secchi, Sheela V. Belur, Javad Ivakpour, Nikolaus Hansen, Alexei Kuntsevich and Franz Kappel, C.T. Kelley, Brecht Donckels, Miroslav Balda, Paul Spencer, Juerg Schwizer, Daniel Buckton, Petr Mikulik, David Gingras, Joachim Vandekerckhove, Yi Cao, Sachin A. Nikumbh, Oliver Bunk, R. G. Abraham, Bruno Luong, J. Ollivier, D. Riley, E. Trautman, F. Esmonde-White';

b = [ 'Matlab/iData version ' vers ' (' date ') by ' auth '.' ];
if nargin > 1
  b = [ b ' Contributions from ' contrib '.' ];
end


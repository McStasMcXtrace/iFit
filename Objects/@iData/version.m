function [b, vers, info] = version(a,long_request)
% v = version(iData): iData class version
%
%   @iData/version: version of the iData class library
%
%   version(iData, 'long')
%     returns the full iData version with all contributors.
%   [str, ver, info] = version(iData)
%     returns as well the short version name, and additional infotmation:
%       memory (total, free, used by Matlab) and number of available CPU's
%
% Version: $Date$ $Author$

% EF 23/09/07 iData implementation

vers = '@IFIT_VERSION@';
date = '@IFIT_DATE@';
auth = '@IFIT_COPYING@';
contrib = 'Eric Ludlam, Felix Morsdorf, Joe Conti, Douglas M. Schwarz, Alexandros Leontitsis, F. Sigworth, Argimiro R. Secchi, Sheela V. Belur, Javad Ivakpour, Nikolaus Hansen, Alexei Kuntsevich and Franz Kappel, C.T. Kelley, Brecht Donckels, Miroslav Balda, Paul Spencer, Juerg Schwizer, Petr Mikulik, David Gingras, Joachim Vandekerckhove, Yi Cao, Oliver Bunk, R. G. Abraham, Bruno Luong, J. Ollivier, D. Riley, E. Trautman, F. Esmonde-White, Y. Altman, Dirk-Jan Kroon, G. Toombes, A. Zheludev, A. Tennant and D. Mc Morrow, Hargreave and Hullah, A. Bouvet/A. Filhol, K. Yamaguchi, W. Falkena, J. Kohlbrecher, J. Rodriguez-Carvajal, W. Baumeister, A. Schmolck, V. Rathod, D. Valevski, J. Almeida, M. Nilsson, J. van Beek, D. Garcia, A. Grinsted, C. Kothe, J. Bialek, D-J Kroon, G. Flandin, M. A. Hopcroft, J. Hokanson, A. Silakov, M. Radin, J. Kantor, C. Pelizzari, J. Yeh, J. Dillon, L. Harriger, G. Romano, Y. Lengwile, W. Falkena, T. Perring, R. Ewings, A. Buts, J. van Duijn, I. Bustinduy, D. Whittaker, S. Toth, D. Alfe, C. Rossant, D. Schwarz, P. Maher, J. Dinale, S. Petit, U. Schwarz, S. Danylenko';

b = [ vers ' iFit/iData (' date ') by ' auth '. $Date$' ];
if nargin > 1
  b = [ b '** Licensed under the EUPL V.1.1 ** Contributions from ' contrib '. Send email to <mailto:ifit-users@mccode.org> to report bugs and requests. More on <http://ifit.mccode.org>.' ];
end

if nargout > 2
  % grab additional information
  info = memoryInfo;
  try; info.NumCores    = feature('NumCores'); end
  try; info.GetOS       = feature('GetOS'); end
  try; info.MatlabPid   = feature('GetPid'); end
  try; info.NumThreads  = feature('NumThreads'); end
  if usejava('jvm') && ~isfield(info,'NumCores')
    r=java.lang.Runtime.getRuntime;
    % mem_avail   = r.freeMemory;
    info.NumCores = r.availableProcessors;
  end
end


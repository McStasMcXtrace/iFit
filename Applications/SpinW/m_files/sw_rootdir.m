function rootdir = sw_rootdir()
% gives the path to the SpinW code
%
% rootdir = SW_ROOTDIR()
%
% See also SW.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

rootdir = mfilename('fullpath');
idx     = strfind(rootdir,filesep);
rootdir = rootdir(1:idx(end-1));

end
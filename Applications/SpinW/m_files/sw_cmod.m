function r = sw_cmod(r, tol)
% modulo one with tolerance
%
% r = SW_CMOD(r, tol)
%
% It calculates modulo one with tolerance, numbers larger than 1-epsilon >
% 1-tol will be converted to -epsilon.
%
% See also MOD.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

if nargin == 0
    help sw_cmod;
    return
end

r = mod(r,1);

r(r > 1-tol) = r(r > 1-tol)-1;

end
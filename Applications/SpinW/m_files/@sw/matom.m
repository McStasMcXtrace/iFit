function mAtomList = matom(obj, varargin)
% generates all magnetic atoms in the unit cell
%
% Same as sw.atom, but only lists the magnetic atoms, which has non-zero
% spin. For speedup, call obj.matom(true), then it reads saved atomic
% positions from obj.
%
% See also SW.ATOM.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

if (nargin>1) && (varargin{1}>0) && ~isempty(obj.matomstore)
    
    mAtomList = obj.matomstore;
else
    
    atomList      = obj.atom;
    
    mAtomList.r   = atomList.r(:,atomList.mag==1);
    mAtomList.idx = atomList.idx(:,atomList.mag==1);
    mAtomList.S   = obj.unit_cell.S(mAtomList.idx);
    
    obj.matomstore = mAtomList;
end

end
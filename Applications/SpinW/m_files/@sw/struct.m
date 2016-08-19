function objS = struct(obj)
% extracts all public properties of sw object into a struct
%
% objS = STRUCT(obj)
%
% See also SW, SW.COPY.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

objS   = struct;
fNames = fieldnames(obj);
for ii = 1:length(fNames)
    objS.(fNames{ii}) = obj.(fNames{ii});
end

end % struct

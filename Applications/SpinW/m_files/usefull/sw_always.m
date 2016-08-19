function out = sw_always(inp)
% converts symbolic logical expressions into logical expressions
%
% out = SW_ALWAYS(inp)
%
% Input:
%
% inp   Any symbolic/logical type matrix.
%
% Output:
%
% out   Logical output with the same dimensions as the input.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

if isa(inp,'sym')
    out = isAlways(inp);
else
    out = logical(inp);
end

end
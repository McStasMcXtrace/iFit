function a = copyobj(a)
% b = copyobj(s) : makes a copy of iData object
%
%   @iData/copyobj function to return a duplicate of data sets.
%   creates a new iData object with same content as 'a', but different Tag/ID and Date.
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=copyobj(a);
%
% Version: $Revision: 1.2 $
% See also iData, iData/uplus, iData/findobj

a = iData_private_history(iData_private_newtag(a), mfilename, a);


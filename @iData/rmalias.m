function s_out = rmalias(a_in,names)
% [s,...] = rmalias(s, AliasName) : removes iData aliases
%
%   @iData/rmalias function to remove iData aliases.
%   The function works also when AliasName is given as a cell string.
%   The command rmalias(iData,'Signal') resets the Signal to the biggest numerical field.
%   The input iData object is updated if no output argument is specified.
%
% input:  s: object or array (iData)
%         AliasName: Name of existing or new alias (char/cellstr)
% output: s: array (iData)
% ex:     rmalias(iData,'Temperature')
%
% Version: $Revision: 1.3 $
% See also iData, iData/getalias, iData/get, iData/set, iData/setalias

% EF 27/07/00 creation
% EF 23/09/07 iData implementation
if nargin == 1
  names='';
end

s_out = setalias(a_in, names);

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),s_out);
end

function this = rmalias(this,names)
% s = rmalias(s, AliasName) : removes iData aliases
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
% Version: $Revision: 1.5 $
% See also iData, iData/getalias, iData/get, iData/set, iData/setalias

% EF 27/07/00 creation
% EF 23/09/07 iData implementation
if nargin == 1
  names='';
end

% handle array of objects
if numel(this) > 1
  for index=1:numel(this)
    this(index) = rmalias(this(index), names);
  end
  return
end

this = setalias(this, names);

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1), this);
end

function b = isfield(a, field)
% b = isfield(s, field) : check existence of field/alias in iData objects
%
%   @iData/isfield function which checks if a name is already defined as a field or alias in the iData object
%
% input:  s: object or array (iData)
%         field: name to check for (string)
% output: b: true when the name is already defined, false otherwise
% ex:     b=isfield(a, 'history');
%
% Version: $Revision: 1.1 $
% See also iData, isfield
b= false;
if any(strcmpi(field, fieldnames(a)))
  b = true;
elseif any(strcmpi(field, getalias(a)))
  b = true;
elseif any(strcmpi(field, {'history','filename','axes','alias','axis'}))
  b = true;
end


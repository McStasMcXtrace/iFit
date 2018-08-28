function b = isfield(a, field)
% b = isfield(s, field) : check existence of field/parameter in model objects
%
%   @iFunc/isfield function which checks if a name is already defined as a Parameter
%     isfield(s) returns true when the field is defined in the object.
%
% input:  s: object or array (iFunc)
%         field: name to check for (string)
% output: b: not empty when the parameter is already defined, false otherwise
% ex:     b=isfield(a, 'Temperature');
%
% Version: $Date$
% See also iFunc, isfield, iFunc/findfield

if nargin == 1, field=[]; end

if isempty(field), b=a.Parameters; return; end

b = any(strcmp(field, a.Parameters));


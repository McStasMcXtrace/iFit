function [b, f] = isfield(a, field)
% b = isfield(s, field) : check existence of field/parameter in model objects
%
%   @iFunc/isfield function which checks if a name is already defined as a Parameter
%     isfield(s, field) returns true when the field is defined in the object.
%
%   [b, match] = isfield(s, field)
%     also returns the field match, which value is: get(s, match)
%
%   A wider/nested search is obtained with: findfield(s, field)
%
% input:  s: object or array (iFunc)
%         field: name to check for (string)
% output: b: true when the parameter is already defined, false otherwise
% ex:     b=isfield(bose, 'Temperature');
%
% Version: $Date$
% See also iFunc, isfield, iFunc/findfield

if nargin == 1, field=[]; end
b=false; f=[]; 

if isempty(field), b=a.Parameters; return; end

for fields={a.Parameters, fieldnames(a)}
  this_fields = fields{1};
  match = strcmp(field, this_fields); % search first for exact match

  if ~any(match)
    match = strcmp(lower(field), lower(strtok(this_fields))); % relaxed search (case insensitive)
  end
  
  if any(match), 
    b = find(match, 1, 'first');
    f = strtok(this_fields{b});
    break
  end
end


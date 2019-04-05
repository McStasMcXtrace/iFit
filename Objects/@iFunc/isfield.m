function [tf, f, pindex] = isfield(a, field)
% tf = isfield(s, field) : check existence of field/parameter in model objects
%
%   @iFunc/isfield function which checks if a name is already defined as a Parameter
%     isfield(s, field) returns true when the field is defined in the object.
%
%   [tf, match,index] = isfield(s, field)
%     also returns the field match, which value is: get(s, match)
%     and the corresponding parameter index.
%
%   A wider/nested search is obtained with: findfield(s, field)
%
% input:  s: object or array (iFunc)
%         field: name to check for (string)
% output: tf: 0 when not found, 1 when is a parameter, 2 when is a class property.
% ex:     tf=isfield(bose, 'Temperature');
%
% Version: $Date$ $Version$ $Author$
% See also iFunc, isfield, iFunc/findfield

if nargin == 1, field=[]; end
tf=false; f=[]; pindex=[];

if isempty(field), tf=a.Parameters; return; end

fields={a.Parameters, fieldnames(a)};
for index=1:numel(fields)
  this_fields = fields{index};
  match = strcmp(field, this_fields); % search first for exact match

  if ~any(match)
    match = strcmp(lower(field), lower(strtok(this_fields))); % relaxed search (case insensitive)
  end
  
  if any(match), 
    tf = index;
    pindex = find(match, 1, 'first');
    f  = strtok(this_fields{tf});
    break
  end
end


function b = isfield(a, field)
% b = isfield(s, field) : check existence of field/alias in iData objects
%
%   @iData/isfield function which checks if a name is already defined as a field or alias in the iData object
%   ifield(s) returns the full list of defined fields and aliases in the object.
%
% input:  s: object or array (iData)
%         field: name to check for (string)
% output: b: true when the name is already defined, false otherwise
% ex:     b=isfield(a, 'history');
%
% Version: $Revision$
% See also iData, isfield

persistent fields

if isempty(fields), fields=fieldnames(iData); end

if nargin == 1, field=[]; end

if numel(a) > 1
  if isempty(field), b = cell(size(a));
  else               b = zeros(size(a)); end
  parfor index=1:numel(a)
    if isempty(field), b{index} = feval(mfilename, a(index), field);
    else               b(index) = feval(mfilename, a(index), field); end
  end
  return
end

if isempty(field)
  b = [ fields ; getalias(a) ; {'history','filename','axes','alias','axis' }' ];
  return
end

b= false;
if any(strcmpi(field, fields))
  b = true;
elseif any(strcmpi(field, getalias(a)))
  b = true;
elseif any(strcmpi(field, {'history','filename','axes','alias','axis'}))
  b = true;
end


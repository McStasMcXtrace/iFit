function [b, link] = fileattrib(a, field, allfields)
% [attribute, link] = fileattrib(s, field) : return a field Attribute
%
%   @iData/fileattrib function which looks for an associated Attribute to a field.
%      Attributes are set from e.g. NetCDF/CDF/NeXus/HDF files.
%      returns []  when no attribute exists
%      returns NaN when the field is already an attribute
%
% input:  s:     object or array (iData)
%         field: Alias/path in the object (string)
% output: attribute: the value of the associated Attribute, or [].
%         link:      the path of the associated Attribute, or [].
% ex:     b=fileattrib(a, 'Signal');
%
% Version: $Revision$
% See also iData, isfield

% handle array of iData input
if numel(a) > 1
  b = cell(1, numel(a)); link=b;
  parfor index=1:numel(a)
    [b{index},link{index}] = feval(mfilename, a(index), field);
  end
  return
end

if nargin < 3
  allfields = findfield(a);
end

if iscellstr(field)  && numel(field) > 1
  b = cell(1, numel(field)); link=b;
  for index=1:numel(field)
    [b{index},link{index}] = feval(mfilename, a, field{index}, allfields);
  end
  return
end

field = char(field);

% replace alias by its link
if any(strcmpi(field, getalias(a)))
  alias = getalias(a, field);
  if ischar(alias)
    field = alias;
  end
end

b = []; link = '';

if nargin == 1
  [status, b] = fileattrib(a.Source);
  link        = a.Source;
  if ~status, b = []; end
else
  
  if any(strcmp(field, allfields))
    [b, link] = iData_getAttribute(a, field);
  end
end


function [b, location] = fileattrib(a, field, varargin)
% FILEATTRIB Get or set field Attributes
%   Looks for an associated Attribute to a field.
%      Attributes are set from e.g. NetCDF/CDF/NeXus/HDF files.
%      returns []  when no attribute exists
%
%   s=FILEATTRIB(s, field, attribute) sets the attribute for given field, when
%   attribute is given as a struct.
%
%   s=FILEATTRIB(s, field, name, value) sets the attribute 'name=value' for given
%   field.
%
% input:  s:     object or array (estruct)
%         field: Alias/path in the object (string)
%         attributes: when given as a structure, sets the attributes for the field.
%         it can also be given as a name/value pair
% output: attribute: the value of the associated Attribute, or [].
%                    or the updated object when storing attributes.
%         link:      the path of the associated Attribute, or [].
% ex:     b=fileattrib(a, 'Signal'); 
%         b=fileattrib(b, 'Signal',struct('long_name','hello world'))
%
% Version: $Date$ $Version$ $Author$
% See also estruct, isfield

% The Attribute property maps the object structure and allows to store any
% information attached to an existing hierarchy.

location = [];
if nargin == 1
  b = get(a, 'Attributes');
  return
end

% handle object arrays
if numel(a) > 1
  b        = {};
  location = {};
  for index=1:numel(a)
    [b{end+1}, location{end+1}] = fileattrib(a(index), field, varargin{:});
  end
  return
end

% handle field array
if iscell(field) && numel(field) > 1
  b        = {};
  location = {};
  for index=1:numel(field)
    [b{end+1}, location{end+1}] = fileattrib(a, field{index}, varargin{:});
  end
  return
end

% first resolve the 'true' path to the field.
if isnumeric(field) && isscalar(field)
  field = getaxis(a, num2str(field));
  if isempty(field), b = []; return; end
end
field = char(field);
try
  f = getalias(a, field);
  if ~isempty(f) && ischar(f) && isfield(a, f)
    if a.verbose > 1
      disp([ mfilename ': field ' field ' is an alias for ' f ]);
    end
    field = f;
  end
end

% we decompose the field into [base].[group].[field] when appropriate.
[base, group, field] = getAttributePath(field);

% Identify existing location of the Attribute for given field.
% Possible attribute locations are listed below.
locations = { [ base group field '.Attributes' ], ...
              [ base group 'Attributes.' field ], ...
              [ 'Attributes.' base group field ], ...
              [ 'Headers.'    base group field ], ...
              [ 'Attributes.'      group field ] };

for loc = locations
  if isfield(a, loc{1}), location = loc{1}; break; end
end

% select default location when does not exist yet.
if isempty(location)
  location = [ 'Attributes.' base group field ];
end
if a.verbose > 1
  disp([ mfilename ': Attribute location ' location ]);
end

try
  b = get(a, location, 'alias');
catch
  b = [];
end
if isempty(b)
  try
    b = get(a, location);
  catch
    if nargin == 2 && a.verbose > 1
      disp([ mfilename ': could not get Attribute ' location ]);
    end
  end
end
  
if nargin == 2
  % get attribute value
  return
else
  % set attribute value
  if nargin == 4 && ischar(varargin{1})
    attribute.(varargin{1}) = varargin{2};
  else
    attribute = varargin{1};
  end
  if isstruct(attribute)
    % we merge with existing
    for ff = fieldnames(attribute)'
      b.(ff{1}) = attribute.(ff{1});
    end
    % and store
    set(a, location, b, 'alias');
    if a.verbose
      disp([ mfilename ': store ' location ]);
    end
  end
end

% ------------------------------------------------------------------------------
function [base, group, lastword] = getAttributePath(field)
% GETATTRIBUTEPATH Cut the entry name into basename, group and dataset
  tokens = textscan(field, '%s', 'Delimiter','./\\');
  tokens = strcat(tokens{1},'.');
  lastword=tokens{end}; lastword(end) = [];

  if numel(tokens) == 1
    base = ''; group = '';
  elseif numel(tokens) == 2
    base = ''; group = tokens{1};
  else
    base = tokens{1}; group = [ tokens{2:(end-1)} ];
  end

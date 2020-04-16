function [b, location] = fileattrib(a, field, varargin)
% FILEATTRIB Get or set field Attributes
%   attribute = FILEATTRIB(s, field) looks for an associated Attribute to a field.
%   Attributes are set from e.g. NetCDF/CDF/NeXus/HDF files. They correspond with
%   additional information attached to object properties.
%   Returns [] when no attribute exists.
%
%   FILEATTRIB(s, field, attribute) sets the attribute for given field.
%   The attribute is often given as a struct. The attribute is returned.
%
%   FILEATTRIB(s, field, name, value) sets the attribute 'name=value' for given
%   field. The attribute is returned.
%
%   [attr, loc] = FILEATTRIB(...) also returns the attribute location.
%
% Example: b=estruct(1:10); fileattrib(b, 'Signal','unit','carrots');...
%          isstruct(fileattrib(b, 'Signal'))
%
% Version: $Date$ $Version$ $Author$
% See also estruct, isfield

% The Attribute property maps the object structure and allows to store any
% information attached to an existing hierarchy.

location = []; b = [];
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
      disp([ mfilename ': DEBUG: field ' field ' is an alias for ' f ]);
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
              [ 'Headers.'    base group field ], ...
              [ 'Attributes.'      group field ], ...
              [ base 'Attributes.' group field ], ...
              [ 'Attributes.' base group field ] };

for loc = unique(locations)
  if isfield(a, loc{1}), location = loc{1}; break; end
end

% select default location when does not exist yet.
if isempty(location)
  location = [ 'Attributes.' base group field ];
end
if a.verbose > 1
  disp([ mfilename ': DEBUG: Attribute location ' location ]);
end

try
  b = get(a, loc{1}, 'alias');
end
if isempty(b)
  try
    b = get(a, loc{1});
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
    if a.verbose > 1
      disp([ mfilename ': DEBUG: store ' location ]);
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

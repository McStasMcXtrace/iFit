function [b, location] = fileattrib(a, field, attribute, value)
% [attribute, link] = fileattrib(s, field) : return a field Attribute
%
%   @estruct/fileattrib function which looks for an associated Attribute to a field.
%      Attributes are set from e.g. NetCDF/CDF/NeXus/HDF files.
%      returns []  when no attribute exists
%      returns NaN when the field is already an attribute
%
%   s=fileattrib(s, field, attribute) sets the attribute for given field, when
%     attribute is given as a struct.
%
%   s=fileattrib(s, field, name, value) sets the attribute 'name=value' for given
%     field.
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

% The Attribute property maps the object structure and allow to store any
% information attached to an existing hierarchy.

location = [];
if nargin == 1
  [~,b] = fileattrib(a.Source);
  return
end

% first resolve the 'true' path to the field.
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
              [ 'Headers.' base group field ]};


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
  if isempty(b)
    b = get(a, location);
  end
catch ME
  b = [];
  if nargin == 2 && a.verbose > 1
    disp([ mfilename ': could not get Attribute ' location ]);
  end
end
  
if nargin == 2
  % get attribute value
  return
else
  % set attribute value
  if nargin == 4 && ischar(attribute)
    sattr.(attribute) = value;
    attribute = sattr;
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
% GETATTRIBUTEPATH the entry name into basename, group and dataset

  % get group and field names
  lastword_index = find(field == '.' | field == '/' | field == '\', 2, 'last'); % get the group and the field name
  if isempty(lastword_index)
    lastword = field; 
    group    = '';
    base     = '';                            % Attributes.<field>.
  elseif isscalar(lastword_index)
    lastword = field((lastword_index+1):end); 
    group    = field(1:lastword_index);
    base     = '';                            % <group>.Attributes.<field>
  else 
    lastword = field( (lastword_index(2)+1):end ); 
    group    = field( (lastword_index(1)+1):lastword_index(2) ); 
    base     = field(1:lastword_index(1));    % <basename>.<group>.Attributes.<field>
  end

  if ~isempty(base)  && base(end)  ~= '.', base = [ base '.' ]; end
  if ~isempty(group) && group(end) ~= '.', group= [ group '.' ]; end

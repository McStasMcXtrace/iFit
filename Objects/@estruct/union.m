function [ai,bi] = union(a, b)
% [ai,bi] = union(a, b) : computes object union area and values
%
%   @estruct/union function to pre-combine data sets
%   This function computes the common axes union between data sets
%   Resulting objects are returned, e.g. for performing further operations
%     ai = union(a) where 'a' is an object array computes union of all elements
%
% input:  a: object or array (estruct)
%         b: object (estruct)
% output: b: object or array (estruct)
% ex:     b=union(a, a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/setaxis, estruct/getaxis, estruct/intersect

if nargin == 2
  bi = union([a b]);
  ai = bi(1);
  bi = bi(2);
  return
end

if numel(a) == 1, ai=a; bi=a; return; end

% first check if all objects have same axes
all_identical_axes=1;
for index=1:ndims(a(1)) % loop on axes
  x = getaxis(a(1), index);
  for obj=2:numel(a)
    if ~isequal(getaxis(a(obj), index), x)
      all_identical_axes=0; break;
    end
  end
  if ~all_identical_axes, break; end
end
% return if using identical axes: no need to interpolate. retain axes and data.
if all_identical_axes, ai=a; bi=[]; return; end

% compute common axes
c_axis = private_caxis(a,'union');

% loop on all estruct to interpolate
ai = a; bi=[];
for index=1:numel(a)
  if ~isempty(a(index))
    ai(index) = interpn(a(index), c_axis);
  end
end


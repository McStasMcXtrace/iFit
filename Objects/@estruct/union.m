function [ai,bi] = union(a, b)
% UNION  Object axes union.
%   [AU,BU]=UNION(A,B) for objects A and B, extends both object axes onto the 
%   union of axes bounds. Resulting objects are returned, e.g. for performing
%   further operations (e.g. combine, plus, ...).
%
%   AU = UNION(A) where 'A' is an object array computes intersection of
%   all elements
%
% Example: a=estruct(peaks); b=copyobj(a); moveaxis(a,1,-5); ...
%   [au,bu]=union(a,b); size(au,1)==54
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


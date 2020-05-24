function [ai,bi] = intersect(a, b)
% INTERSECT Object axes intersection.
%   [AI,BI]=INTERSECT(A,B) for objects A and B, returns the subsets common to
%   the two objects axes. Resulting objects are returned, e.g. for performing
%   further operations (e.g. combine, plus, ...).
%
%   AI = INTERSECT(A) where 'A' is an object array computes intersection of
%   all elements
%
% Example: a=iData(peaks); b=copyobj(a); moveaxis(a,1,-5); ...
%   [ai,bi]=intersect(a,b); size(ai,1)==44
% Version: $Date$ $Version$ $Author$
% See also iData, iData/setaxis, iData/getaxis, iData/interp, iData/union

if nargin == 2
  bi = intersect([a b]);
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
c_axis = private_caxis(a,'intersection');

% loop on all iData to interpolate
ai = zeros(iData, numel(a),1); bi=[];
for index=1:numel(a)
  ai(index) = interpn(a(index), c_axis(1:ndims(a(index))));
end


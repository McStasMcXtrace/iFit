function b = meshgrid(a, varargin)
% s = meshgrid(a) : transforms an iData object so that its axes are grids
%
%   @iData/meshgrid function to transform iData object axes so that they are
%     on a regular grid, as obtained from ndgrid. When the initial axes are not
%     perpendicular/regular, the object Signal is interpolated on the new grid.
%   meshgrid(a, 'vector' ...) forces all axes as vectors
%   meshgrid(a, dims, ...)    specifies the size of the histogram
%   meshgrid(a, 'fill' ...)   fills empty histogram bins when converting from an event list
%
%   A meshgrid histogram can be converted into an event list with the 'event' method.
%
% input:  a: object or array (iData)
%         method: 'linear','cubic','spline','nearest'
%                 'vector' to get only vector axes
% output: s: object (iData)
% ex:     c=meshgrid(a); c=meshgrid(a, 'vector linear')
%         c=meshgrid(a, 100, 'fill')
%
% Version: $Date$
% See also iData, iData/interp, iData/hist, iData/event

% handle input iData arrays
if numel(a) > 1
  b = zeros(iData, numel(a), 1);
  parfor index=1:numel(a)
    b(index) = meshgrid(a(index),varargin{:});
  end
  b = reshape(b, size(a));
  return
end

method = '';
n_dims = [];
if nargin > 1
  for index=1:numel(varargin)
    if ischar(varargin{index})
      method = [ method ' ' varargin{index} ];
    elseif isnumeric(varargin{index})
      n_dims = varargin{index};
    end
  end
end

s_dims = size(a); % Signal/object dimensions

if isscalar(n_dims)
  n_dims = n_dims*ones(1,ndims(a));
end

% check axes dimensions
for index=1:ndims(a)
  v = getaxis(a, index);

  % compute the initial axis length
  if isvector(v), a_len = numel(v);
  else            a_len = size( v, index);
  end
  if isvector(a) >= 2 && a_len > prod(size(a))^(1/ndims(a))*2 % event data set
    a_len = ceil(prod(size(a))^(1/ndims(a))*2);
  end
  if a_len == 1, a_len = 2; end
  if numel(n_dims) < index
    n_dims(index) = a_len;
  end
end
clear v
if ndims(a) == 1
  n_dims = [ max(n_dims) 1 ];
end

if isvector(a) > 2
  b = hist(a, n_dims, method);
  return;
end

% create a regular grid
[f_axes, changed] = iData_meshgrid(a, n_dims, method);
method            = strtrim(strrep(method, 'vector', ''));

b = copyobj(a);
% interpolate if needed
if changed
  b = interp(b, f_axes, method);
else
  % transfer grid axes to object
  for index=1:length(f_axes)
    b = setaxis(b, index, f_axes{index});
  end
end

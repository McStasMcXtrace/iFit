function b = meshgrid(a, varargin)
% MESHGRID Object rebinned onto grid axes.
%   B = MESHGRID(A) transforms object axes so that they are on a regular grid, 
%   as obtained from ndgrid. When the initial axes are not perpendicular/regular,
%   the object Signal is interpolated on the new grid.
%   A meshgrid histogram can be converted back into an event list with the EVENT 
%   method.
%
%   B = MESHGRID(A, [M N ...]) specifies the size of the rebinned histogram. When
%   given as a single scalar value M, the same dimension is used for all axes.
%
%   B = MESHGRID(..., 'vector') requests all axes to be set as vectors.
%
%   B = MESHGRID(..., 'grid') requests all axes as grid with same size as Signal
%
%   B = MESHGRID(..., 'fill') fills empty histogram bins (e.g. when converting 
%   from an event list) will FILL so that the final object has no NaN's and empty
%   bins.
%
%   B = MESHGRID(..., METHOD ...) uses specified METHOD for interpolation as one of
%   'linear' (default), 'spline', 'cubic', or 'nearest'.
%
% Example: a=iData(peaks); b=meshgrid(a, 101); size(b,1)==101
% Version: $Date$ $Version$ $Author$
% See also iData, iData/interp, iData/hist, iData/event, iData/fill

% handle input iData arrays
if numel(a) > 1
  b = zeros(iData, numel(a), 1);
  for index=1:numel(a)
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
myisvector = @(c)length(c) == numel(c);
for index=1:ndims(a)
  v = getaxis(a, index);

  % compute the initial axis length
  if myisvector(v), a_len = numel(v);
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
[f_axes, changed] = private_meshgrid(a, n_dims, method);
fillme            = strfind(method, 'fill');
method            = strtrim(strrep(method, 'vector', ''));

b = copyobj(a);
% interpolate if needed
if changed
  b = interpn(b, f_axes, method);
else
  % transfer grid axes to object
  for index=1:length(f_axes)
    b = setaxis(b, index, f_axes{index});
  end
end

if ~isempty(fillme)
  b = fill(b);
end

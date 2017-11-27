function b = event(a, dim, signal, err, monitor)
% c = event(a) : build an event list 
%
%   @iData/event function to transform an object into an event list.
%
%     An nD event list object is a list of points (e.g. one per row) for
%       which the points coordinates are given as the n first columns
%       and subsequent columns define a Signal.
%
%     event(a) will attempt to guess the dimensionality (up to 3D)
%     event(a, dim) specifies the object dimensionality. The first dim
%       columns are used for the coordinates. 'dim' is an integer.
%     event(a, dim, signal) also specifies which column has to be used as
%       Signal. When given as a vector, the norm of the specified columns is used.
%     event(a, dim, signal, error, monitor) also specifies which column
%       has to be used as Error and Monitor.
%
%     An event list object can be converted back into a mesh-type histogram
%     with the meshgrid or hist method, e.g.
%       b=event(a); c=meshgrid(b);
%     then a and c are close.
%
% input:  a:   object or array (iData)
%         m,n,p...: dimensions (integers)
% output: c: object or array (iData)
% ex:     a=iData(rand(50,4)); b=event(a);
%
% Version: $Date$
% See also iData, iData/squeeze, iData/size, iData/permute, iData/resize, iData/meshgrid, iData/hist

% first get dimensions from varargin
if nargin < 2, dim=[]; end
if nargin < 3, signal=[]; end
if nargin < 4, err=[]; end
if nargin < 5, monitor=[]; end

% handle iData array
if numel(a) > 1
  for index=1:numel(a)
    a(index) = feval(mfilename, a(index), dim, signal, err, monitor);
  end
  if nargout == 0 && ~isempty(inputname(1))
    assignin('caller',inputname(1),a);
  end
  return
end


if ndims(a) > 2 || size(a, 2) > 20
  % case: we have an histogram (not a single vector)
  
  % create an object with meshgrid type axes plus Signal
  
  % we test if the axes are all grids with same size as the signal.
  % In this case we just need to reshape all of the stuf as long vectors
  sig = getaxis(a, 0); sig = numel(sig);
  flag = false;
  for index=1:ndims(a)
    x = getaxis(a, index);
    if numel(x) ~= sig
      flag = true;  % need to meshgrid
      break
    end
  end
  
  if flag % need to meshgrid
    b = meshgrid(a,'grid');
  else
    b = copyobj(a);
  end
  for index=1:ndims(a)
    x = getaxis(a, index);
    xl= a.Alias.Axis{index};
    if ~ischar(xl), xl = sprintf('Axis_%i', index); end
    b = setalias(b, [ xl '_event' ], x(:));
    setaxis(b, index, [ xl '_event' ]);
  end
  % handle Signal, Error and Monitor
  x = get(b, 'Signal'); 
  b = set(b, 'Signal', x(:));
  x = get(b, 'Error'); 
  b = set(b, 'Error', x(:));
  x = get(b, 'Monitor'); 
  b = set(b, 'Monitor', x(:));
  
elseif ndims(a) == 2 && size(a, 2) <= 20 % e.g [x,y,z, bx,by,bz]
  % case: it really looks like a 2D column-type data set

  % try to guess what to do, if the dimensionality is not given
  if isempty(dim)
    dim = min(3, size(a,2)-1);
  end
  if isempty(signal), signal=dim+1; end

  data = get(a, 'Signal');
  b = copyobj(a);
  if ~isempty(signal) && isnumeric(signal)
    s = data(:, signal)';
    if all(size(s) > 1), s = sqrt(sum(s .* s)); end
    b = set(b, 'Signal', s(:));
  end
  if isnumeric(err) && numel(err) == 1
    b = set(b, 'Error', data(:, err));
  end
  if isnumeric(monitor) && numel(monitor) == 1
    b = set(b, 'Monitor', data(:, monitor));
  end
  % set axes
  for index=1:dim
    b = setaxis(b, index, data(:,index));
  end
end

% set return value
if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),b);
end

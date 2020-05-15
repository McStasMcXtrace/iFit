function b = interpn(a, varargin)
% INTERPN N-D object interpolation (table lookup).
%    VI = INTERPN(V, Y1,Y2,Y3,...) interpolates object V onto new axes Y1,Y2,Y3, 
%    etc to find the underlying object VI values.
%    INTERPN works for all object dimensionalities. The multidimensional interpolation
%    is based on a Delauney tessellation using the Computational Geometry 
%    Algorithms Library, CGAL. Extrapolated data is set to NaN for the Signal, 
%    Error and Monitor. Also, in some cases, the triangulation interpolant creates
%    fake 'flat' area, especially when the axes area is concave. We then recommend
%    you try the HIST method.
%    For Event data sets, we recommend to use the HIST method which is much faster.
%
%    VI = INTERPN(V, {Y1,Y2,Y3,...}) is similar to the previous syntax.
%
%    VI = INTERPN(V, A) where A is an estruct object interpolates V onto A axes,
%    i.e. Y1=A{1}, Y2=A{2}, ... with the usual estruct axis notation.
%
%    VI = INTERPN(...,METHOD) specifies alternate methods.  The default
%    is linear interpolation.  Available methods are:
% 
%      'nearest' - nearest neighbor interpolation
%      'linear'  - linear interpolation (default)
%      'spline'  - spline interpolation
%      'cubic'   - cubic interpolation as long as the data is uniformly
%                  spaced, otherwise the same as 'spline'
%
%    VI = INTERPN(..., 'grid') uses meshgrid/ndgrid to determine new axes as arrays.
%
%    VI = INTERPN(..., 'vector') requests all axes to be set as vectors.
%
% Example: a=estruct(peaks); b=interpn(a, 'grid'); isequal(a,b)
% Version: $Date$ $Version$ $Author$
% See also estruct, interp1, interpn, ndgrid, estruct/setaxis, estruct/getaxis,
%          estruct/hist, estruct/resize, estruct/reshape, estruct/fill

% private_interp and private_meshgrid are in private

% handle input arrays
if numel(a) > 1
  b = zeros(estruct, numel(a), 1);
  for index=1:numel(a)
    b(index) = interpn(a(index), varargin{:});
  end
  b = reshape(b, size(a));
  return
end

% build new estruct object to hold the result
if isempty(a), axescheck(a); end % make sure we have something to interpolate (check first)

% object check
if isempty(a), b=a; return; end

% default axes/parameters
i_axes = cell(1,ndims(a)); i_labels=i_axes;
for index=1:ndims(a)
  [i_axes{index}, i_labels{index}] = getaxis(a, index);  % loads object axes, or 1:end if not defined
end
for index=ndims(a):numel(a.Axes)
  [~, i_labels{index}] = getaxis(a, index);  % additional inactive axes labels (used to create new axes)
end
method='linear';

% interpolation axes
f_axes           = i_axes;
requires_meshgrid= 0;

% parse varargin to overload defaults and set manually the axes ----------------
axis_arg_index   = 0;
for index=1:length(varargin)
  c = varargin{index};
  if ischar(c) && ~isempty(strfind(c,'grid'))
    requires_meshgrid=1;
  elseif ischar(c) && ~isempty(c)         % method (char)
    method = c;
  elseif isa(varargin{index}, 'estruct')  % object axes
    if length(c) > 1
      if a.verbose
        warning(['%s: Can not interpolate onto all axes of %i-th input argument which is an array of [%i] elements.\n' ...
        '\tUsing first element only.'], ...
        mfilename, index, numel(c));
      end
      c = c(1);
    end
    for j1 = 1:ndims(c)
      axis_arg_index = axis_arg_index+1;
      [f_axes{axis_arg_index}, lab] = getaxis(c, j1);
      % update labels when missing with that from 2nd object
      if ~isempty(lab) && axis_arg_index < length(i_labels) && isempty(i_labels{axis_arg_index})
        i_labels{axis_arg_index} = lab;
      end
    end
  elseif isnumeric(c) && length(c) ~= 1  % vector/matrix axes
    axis_arg_index = axis_arg_index+1;
    if ~isempty(c), f_axes{axis_arg_index} = c; end
  elseif iscell(c)                      % cell(vector/matrix) axes
    for j1 = 1:length(c(:))
      axis_arg_index = axis_arg_index+1;
      if ~isempty(c{j1}), f_axes{axis_arg_index} = c{j1}; end
    end
  elseif ~isempty(c) && a.verbose
    warning([ mfilename ': Input argument ' num2str(index) ' of class ' class(c) ' size [' num2str(size(c)) '] is not supported. Ignoring.']);
  end
  clear c
end % input arguments parsing

% check for method to be valid
if ~any(strcmp(method, {'linear','cubic','spline','nearest','v4'}))
  if a.verbose
    warning([ mfilename ': Interpolation method "' method '" is not supported. Use: linear, cubic, spline, nearest, v4. Defaulting to linear.']);
  end
  method = 'linear';
end

% check/determine output axes for interpolation --------------------------------

% test axes and decide to call meshgrid if necessary

if isvector(a) >=2 % event data set: redirect to hist method (accumarray)
  f_axes = private_meshgrid(f_axes, s_dims, 'vector'); % private function
  b = hist(a, f_axes{:});
  return
end

% check final axes
i_dims = size(a); % Signal/object dimensions (initial)
f_dims = i_dims;  % Signal/object dimensions (final)
myisvector = @(c)length(c) == numel(c);

% do we have mixed vector and matrix i_axes initial axes ?
has_vector = 0;
has_matrix = 0;
for index=1:ndims(a)
  v = f_axes{index};
  if isempty(v), v= i_axes{index}; end % no axis specified, use the initial one

  % compute the initial axis length
  if myisvector(v), a_len = numel(v); has_vector=1;
  else            a_len = size( v, index); has_matrix=1;
  end
  if isvector(a) >= 2 && a_len > prod(size(a))^(1/ndims(a))*2 % event data set
    a_len = prod(size(a))^(1/ndims(a))*2;
  end
  if a_len == 1, a_len = 2; end
  f_dims(index) = a_len;
end
if length(f_axes) > 1 && ((has_vector && has_matrix) || requires_meshgrid)
  i_axes = private_meshgrid(i_axes, i_dims, method); % private function
end

% do we have mixed vector and matrix f_axes final axes ?
has_vector = 0;
has_matrix = 0;
for index=1:ndims(a)
  v = f_axes{index};
  if isempty(v), v= i_axes{index}; end % no axis specified, use the initial one

  % compute the initial axis length
  if myisvector(v), a_len = numel(v); has_vector=1;
  else            a_len = size( v, index); has_matrix=1;
  end
  if isvector(a) >= 2 && a_len > prod(size(a))^(1/ndims(a))*2 % event data set
    a_len = prod(size(a))^(1/ndims(a))*2;
  end
  if a_len == 1, a_len = 2; end
  f_dims(index) = a_len;
end
if has_vector && has_matrix, requires_meshgrid=1; end

% check if interpolation is indeed required ------------------------------------
% do we need to recompute the final axes ?
if length(f_axes) > 1 && requires_meshgrid
  f_axes = private_meshgrid(f_axes, f_dims, method); % private function
end

% test if interpolation axes have changed w.r.t input object (for possible quick exit)
has_changed = 0;

for index=1:ndims(a)
  this_i = i_axes{index}; if myisvector(this_i), this_i=this_i(:); end
  this_f = f_axes{index}; if myisvector(this_f), this_f=this_f(:); end
  if ~isequal(this_i, this_f)
    % length changed ?
    if length(this_i) ~= length(this_f)
      % not same length
      has_changed=1;
    elseif prod(size(this_i)) ~= prod(size(this_f)) % nb of elements has changed, including matrix axes ?
      has_changed=2;
    elseif all(abs(this_i(:) - this_f(:)) > 1e-4*abs(this_i(:) + this_f(:))/2)
      % or axis variation bigger than 0.01 percent anywhere
      has_changed=3;
    end
  end
  clear this_i this_f
end

% get Signal, error and monitor.
i_signal = subsref(a,struct('type','.','subs','Signal'));

% quick exit check based on the Signal
if any(isnan(i_signal(:))), has_changed=1; end
if ~has_changed && ~requires_meshgrid
  b = a;
  return;
end

if isvector(i_signal) % force axes to be vectors
  for index=1:length(i_axes)
    x=i_axes{index}; x=x(:); i_axes{index}=x;
  end
  clear x
end

% prepare interpolation Signal, Error, Monitor ---------------------------------
i_class    = class(i_signal); i_signal = double(i_signal);

i_error = getalias(a, 'Error');
if ~isempty(i_error),
  % check if Error is sqrt(Signal) or a constant
  if strcmp(i_error, 'sqrt(this.Signal)')
    i_error=[];
  elseif isnumeric(i_error) && isscalar(i_error) == 1
    % keep that as a constant value
  else
    % else get the value
    i_error  = subsref(a,struct('type','.','subs','Error'));
  end
  i_error    = double(i_error);
end

i_monitor = getalias(a, 'Monitor');
if ~isempty(i_monitor),
  % check if Monitor is 1 or a constant
  if isnumeric(i_monitor) && isscalar(i_monitor) == 1
    % keep that as a constant value
  else
    % else get the value
    i_monitor  =subsref(a,struct('type','.','subs','Monitor'));
  end
  i_monitor    = double(i_monitor);
end

% check f_axes vector orientation
for index=1:ndims(a)
  i_axes{index} = double(i_axes{index});
  f_axes{index} = double(f_axes{index});
  if myisvector(f_axes{index})
    % orient the vector along the dimension
    n = ones(1,ndims(a));
    n(index) = numel(f_axes{index});
    if length(n) == 1, n=[ n 1]; end
    f_axes{index}=reshape(f_axes{index},n);
  end
end

% make sure input axes are monotonic. output axes should be OK ------------
i_nonmonotonic=0;
for index=1:ndims(a)
  dif = diff(i_axes{index},1,index);
  if all(dif(:) < 0) || all(dif(:) >= 0), continue; end
  i_nonmonotonic=index; break;
end

if i_nonmonotonic && length(i_axes) > 1
  % transform the initial data into individual points, then interpolate on
  % a regular grid
  i_axes_new  = private_meshgrid(i_axes, size(a));
  % must make sure initial axes given as vector have same length as signal
  flag_ndgrid_needed = 0;
  for index=1:length(i_axes)
    x=i_axes{index}; x=x(:);
    % make sure length(x) == numel(i_signal) for interpolation to work
    if length(x) < numel(i_signal)
      if isempty(length(x) == size(i_signal))
        error('%s: The axis rank %d of length=%d does not match the object %s Signal dimensions %s\n', ...
          mfilename, index, length(x), a.Tag, mat2str(size(i_signal)) ...
        );
      else
        flag_ndgrid_needed = 1;
      end
    end
    i_axes{index} = x;  % now a vector...
  end
  % signal is a grid but axes are vectors, axes should also be...
  if flag_ndgrid_needed
    [i_axes{:}] = ndgrid(i_axes{:});
    for index=1:length(i_axes)
        i_axes{index} = i_axes{index}(:);
    end
  end

  % interpolate initial data on monotonic initial axes
  i_signal    = private_interp(i_axes, i_signal(:),  i_axes_new, method);

  if isnumeric(i_error) && length(i_error) > 1,
    i_error   = private_interp(i_axes, i_error(:),   i_axes_new, method);
  end
  if isnumeric(i_monitor) && length(i_monitor) > 1,
    i_monitor = private_interp(i_axes, i_monitor(:), i_axes_new, method);
  end
  i_axes = i_axes_new;
  clear i_axes_new
end

% last test to check if axes have changed ---------------------------------
has_changed = 0;
for index=1:ndims(a)    % change to double before interpolation
  i_axes{index}=double(i_axes{index});
  f_axes{index}=double(f_axes{index});
end
for index=1:ndims(a)
  x = i_axes{index}; x=x(:)';
  if ~isequal(i_axes{index}, f_axes{index})
    has_changed = 1;
    break
  end
end
if ~has_changed,
  b=a; return;
end

% interpolation takes place here ------------------------------------------
[f_signal, meth] = private_interp(i_axes, i_signal, f_axes, method);
if isscalar(f_signal), b=a; return; end % single scalar value

if isnumeric(i_error) && length(i_error) > 1,
     f_error = private_interp(i_axes, i_error,  f_axes, method);
else f_error = i_error; end
clear i_error
if isnumeric(i_monitor) && length(i_monitor) > 1,
     f_monitor = private_interp(i_axes, i_monitor,f_axes, method);
else f_monitor = i_monitor; end
clear i_monitor
clear i_axes

% get back to original Signal class
if ~strcmp(i_class, 'double')
  f_signal = feval(i_class, f_signal);
  f_error  = feval(i_class, f_error);
  f_monitor= feval(i_class, f_monitor);
end

if isvector(i_signal) && size(i_signal,1)==1
    f_signal = transpose(f_signal);
    f_error  = transpose(f_error);
    f_monitor= transpose(f_monitor);
end
clear i_signal

% transfer Data and Axes --------------------------------------------------
b = zeros(a);  % new object derived from 'a'
set(b, 'Signal' , f_signal);  clear f_signal
set(b, 'Error',   f_error);   clear f_error
set(b, 'Monitor', f_monitor); clear f_monitor
for index=1:length(f_axes)
  set(b, [ 'axis' num2str(index) ], f_axes{index});
  setaxis(b, index, [ 'axis' num2str(index) ]);
end
setalias(b, 'Interpolation', meth);
label(   b, 'Interpolation', 'Interpolation method used');
history(b, mfilename, a, varargin{:});


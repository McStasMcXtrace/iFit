function b = interp(a, varargin)
% [b...] = interp(s, ...) : interpolate iData object
%
%   @iData/interp function to interpolate data sets
%   This function computes the values of the object 's' interpolated
%   on a new axis grid, which may be specified from an other object, as independent axes,
%   or as a rebinning of the original axes.
%     b=interp(s)    rebin/check 's' on a regular grid.
%     b=interp(s, d) where 'd' is an iData object computes 's' on the 'd' axes.
%     b=interp(s, X1,X2, ... Xn) where 'X1...Xn' are vectors or matrices as obtained 
%                    from ndgrid computes 's' on these axes.
%     b=interp(s, {X1,X2, ... Xn}) is similar to the previous syntax
%     b=interp(s, ..., ntimes) where 'ntimes' is an integer computes new axes for 
%                    interpolation by sub-dividing the original axes ntimes.
%     b=interp(s, ..., 'method') uses specified method for interpolation as one of
%                    linear (default), spline, cubic, or nearest
%     b=interp(s, ..., 'grid') uses meshgrid/ndgrid to determine new axes as arrays
%     b=interp(s, ..., 'nan')  also interpolates on nan values
%   Extrapolated data is set to 0 for the Signal, Error and Monitor.
%
% input:  s: object or array (iData)
%         d: single object from which interpolation axes are extracted (iData)
%            or a cell containing axes d={X1,X2, ... Xn}               (cell)
%         X1...Xn: vectors or matrices specifying axis for 
%            dimensions 1 to ndims(s) (double vector/matrix)
%         ntimes: original axis sub-division (integer)
% output: b: object or array (iData)
% ex:     b=interp(a, 'grid','nan');
%
% Version: $Revision: 1.28 $
% See also iData, interp1, interpn, ndgrid, iData/setaxis, iData/getaxis

% input: option: linear, spline, cubic, nearest
% axes are defined as rank of matrix dimensions
% plot function is plot(y,x,Signal)
% rand(10,20) 10 rows, 20 columns
% pcolor/surf with view(2) shows x=1:20, y=1:10

% handle input iData arrays
if length(a) > 1
  b = a;
  for index=1:length(a(:))
    b(index) = interp(a(index), varargin{:});
  end
  return
end

% build new iData object to hold the result
cmd=a.Command;
b = copyobj(a);

% object check
if ndims(a) == 0
  iData_private_warning(mfilename,['Object ' inputname(1) ' ' a.Tag ' is empty. Nothing to interpolate.']);
  return
end
% removes warnings during interp
try
  warn.set = warning('off','iData:setaxis');
  warn.get = warning('off','iData:getaxis');
catch
  warn = warning('off');
end
% default axes/parameters
i_axes = cell(1,ndims(a)); i_labels=i_axes;
for index=1:ndims(a)
  [i_axes{index}, i_labels{index}] = getaxis(a, index);  % loads object axes, or 1:end if not defined 
end
for index=ndims(a):length(a.Alias.Axis)
  [dummy, i_labels{index}] = getaxis(a, index);  % additional inactive axes labels (used to create new axes)
end
method='linear';
ntimes=0;

% interpolation axes
f_axes           = i_axes;
requires_meshgrid= 0; has_grid_arg=0; has_nan_arg=0;

% parse varargin to overload defaults and set manually the axes
axis_arg_index   = 0;
for index=1:length(varargin)
  c = varargin{index};
  if ischar(c) & ~isempty(strfind(c,'grid')) 
    requires_meshgrid=1; has_grid_arg=1;
  elseif ischar(c) & ~isempty(strfind(c,'nan')) 
    has_nan_arg=1;
  elseif ischar(c)                      % method (char)
    method = c;
  elseif isa(varargin{index}, 'iData')  % set interpolation axes: get axis from other iData object
    if length(c) > 1
      iData_private_warning(mfilename,['Can not interpolate onto all axes of input argument ' num2str(index) ' which is an array of ' num2str(numel(c)) ' elements. Using first element only.']);
      c = c(1);
    end
    for j1 = 1:ndims(c)
      axis_arg_index = axis_arg_index+1;
      [f_axes{axis_arg_index}, lab] = getaxis(c, j1);
      if ~isempty(lab) && axis_arg_index < length(i_labels) && isempty(i_labels{axis_arg_index})
        i_labels{axis_arg_index} = lab;
      end
    end
  elseif isnumeric(c) & length(c) ~= 1   % set interpolation axes: vector/matrix
    axis_arg_index = axis_arg_index+1;
    if ~isempty(c), f_axes{axis_arg_index} = c; end
  elseif isnumeric(c) & length(c) == 1  % ntimes rebinning
    ntimes=c;
  elseif iscell(c)                      %set interpolation axes: cell(vector/matrix)
    for j1 = 1:length(c(:))
      axis_arg_index = axis_arg_index+1;
      if ~isempty(c{j1}), f_axes{axis_arg_index} = c{j1}; end
    end
  elseif ~isempty(c)
    iData_private_warning(mfilename,['Input argument ' num2str(index) ' of class ' class(c) ' size [' num2str(size(c)) '] is not supported. Ignoring.']);
  end
end

% check for method to be valid
if isempty(any(strcmp(method, {'linear','cubic','spline','nearest'})))
  iData_private_error(mfilename,['Interpolation method ' method ' is not supported. Use: linear, cubic, spline, nearest.']);
end

% test axes and decide to call meshgrid if necessary
is_grid=0;
if isvector(a) >= 2, requires_meshgrid=1; end % plot3/event style
if ndims(a) > 1
  for index=1:ndims(a)
    % test for the target axes in case they are given as scalars (axes spacing)
    if isscalar(f_axes{index}) && ~isscalar(i_axes{index})
      x = i_axes{index}; x=unique(x); % also makes it a vector
      f_axes{index} = min(x):f_axes{index}:max(x);
    end
    if isvector(f_axes{index}) % vectors should be oriented the right way
      d=ones(1, ndims(a));
      d(index) = length(f_axes{index});
      f_axes{index} = reshape(f_axes{index}, d);
    end
    % this axis is a vector, but others are grids: require meshgrid.
    if any(size(f_axes{index}) == 1)       & is_grid, requires_meshgrid=1; end 
    try
      % this axis is a grid, others should also be...
      if all(size(f_axes{index}) == size(a)),  is_grid=is_grid+1; end 
    end
  end
end

% trigger regular axis check/rebin when interp(a) called
if nargin == 1 & ~is_grid, ntimes=1; end

if ntimes ~= 0
  % rebin iData object using the smallest axes steps for new axes
  for index=1:ndims(a)
    x = i_axes{index}; x=unique(x); % also makes it a vector
    a_step = diff(x);
    a_step = a_step(find(a_step));
    a_step = min([mean(abs(a_step)) median(abs(a_step)) ]);  % smallest non-zero axis step
    if (a_step < 0), a_step = (a_max - a_min)/length(x); end
    a_min  = min(x);
    a_max  = max(x);
    a_len  = (a_max - a_min)/a_step;
    if ntimes > 0 & ntimes < 10
      a_len = a_len*ntimes;
    else
      a_len  = min(a_len, length(x)*10); % can not reduce or expand more 
                                           % than 10 times each axis
    end
    f_axes{index} = linspace(a_min,a_max,a_len+1);
  end
end

% some of the axes are non consistent or not grid: re-bin onto vector axes
if ~is_grid | (requires_meshgrid & is_grid ~= ndims(a))
  % first make axes unique as vectors (sorted)
  for index=1:ndims(a)
    % make the axis as a vector on each dimension
    s = size(a);
    n = ones(1, ndims(a)); if ndims(a) == 1, n = [n 1 ]; end
    if length(find(size(a) > 1)) == 1, n(index) = max(size(a)); % plot3 like
    else n(index) = s(index); end
    v = f_axes{index}; 
    v = unique(v(:));
    f_axes{index} = v;  % vector
  end
end

% test if interpolation axes have changed w.r.t input object
has_changed = 0;
for index=1:ndims(a)  
  this_a = i_axes{index}; if isvector(this_a), this_a=this_a(:); end
  this_i = f_axes{index}; if isvector(this_i), this_i=this_i(:); end
  if ~isequal(this_a, this_i)
    % length changed ?
    if length(this_a) ~= length(this_i)
      % not same length
      has_changed=1; 
    elseif prod(size(this_a)) ~= prod(size(this_i)) % nb of elements has changed, including matrix axes ?
      has_changed=1; 
    elseif all(abs(this_a(:) - this_i(:)) > 1e-4*abs(this_a(:) + this_i(:))/2)
      % or axis variation bigger than 0.01 percent anywhere
      has_changed=1;
    end
  end
end

% get Signal, error and monitor.
i_signal   = get(a,'Signal');
if any(isnan(i_signal(:))) && has_nan_arg, has_changed=1; end
if ~has_changed & (~requires_meshgrid | is_grid), return; end

i_class    = class(i_signal); i_signal = double(i_signal);

i_error = getalias(a, 'Error');
if ~isempty(i_error),   
  % check if Error is sqrt(Signal) or a constant
  if strcmp(i_error, 'sqrt(this.Signal)')
    i_error=[];
  elseif isnumeric(i_error) && isscalar(i_error) == 1
    % keep that as  a constant  value
  else
    % else get the value
    i_error  =get(a,'Error');
  end
  i_error    = double(i_error);
end
  
i_monitor = getalias(a, 'Monitor');
if ~isempty(i_monitor),   
  % check if Monitor is 1 or a constant
  if isnumeric(i_monitor) && isscalar(i_monitor) == 1
    % keep that as a constant  value
  else
    % else get the value
    i_monitor  =get(a,'Monitor');
  end
  i_monitor    = double(i_monitor);
end

% do we need to compute a new grid from axis ?
if ndims(a) > 1 && (length(i_signal) ~= numel(i_signal))
  requires_meshgrid = 1;
end

% compute a new grid before interpolating
if requires_meshgrid
  if length(f_axes) == 1 || (ndims(a) == 1 && isvector(a) == 1)
    % nothing to do as we have only one axis, no grid
  else
    % call ndgrid
    [f_axes{1:ndims(a)}] = ndgrid(f_axes{:});
    if ~has_grid_arg % reshape axes as vectors (but not for 'grid')
      for index=1:ndims(a) 
        f_axes{index} = unique(f_axes{index});
        n = ones(1,ndims(a));
        n(index) = length(f_axes{index});
        if length(n) == 1, n=[ n 1]; end
        f_axes{index}=reshape(f_axes{index},n);
      end
    end
  end
end
% handle 'nan' option: remove the values from the axes/signal, and make it nD vectors
if has_nan_arg
  i_signal = i_signal(:);
  index = find(~isnan(i_signal));
  i_signal = i_signal(index);
  if length(i_axes) == 1 || (ndims(a) == 1 && isvector(a) == 1)
    % nothing to do as we have only one axis, no grid
  else
    % call ndgrid
    [i_axes{1:ndims(a)}] = ndgrid(i_axes{:});
    for i=1:ndims(a)
      x = i_axes{i}; x=x(:); 
      x = x(index); 
      i_axes{i} = x;
    end
  end
end

% make sure input axes are monotonic. output axes should be OK.
i_nonmonotonic=0;
for index=1:ndims(a)
  if any(diff(i_axes{index},1,index) <= 0)
    i_nonmonotonic=1; break;
  end
end

if i_nonmonotonic && ~has_nan_arg
  for index=1:ndims(a)  % apply unique on axes and reorder signal
    i_idx{index}=1:size(a, index);
    [i_axes{index}, i_idx{index}] = unique(i_axes{index});
    if length(i_idx{index}) ~= size(a,index)
      for j=1:ndims(a), 
        if j ~= index, f_idx{j}=':';
        else           f_idx{j}=i_idx{index}; end
      end
      i_signal =i_signal(f_idx{:});
      if isnumeric(i_error) && length(i_error) > 1, 
          try   i_error  =i_error(f_idx{:});
          catch i_error=[]; end
      end
      if isnumeric(i_error) && length(i_monitor) > 1, 
          try   i_monitor=i_monitor(f_idx{:});
          catch i_monitor=[]; end
      end
    end
  end
end

% last test to check if axes have changed
has_changed = 0;
for index=1:ndims(a)
  if ~isequal(i_axes{index}, f_axes{index})
    has_changed = 1;
    break
  end
end
if ~has_changed, return; end

% interpolation takes place here
f_signal = iData_interp(i_axes, i_signal, f_axes, method);
if isnumeric(i_error) && length(i_error) > 1, 
  f_error  = iData_interp(i_axes, i_error,  f_axes, method); 
else f_error = i_error; end
if isnumeric(i_monitor) && length(i_monitor) > 1, 
  f_monitor= iData_interp(i_axes, i_monitor,f_axes, method);
else f_monitor = i_monitor; end

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

% transfer Data and Axes
Data = a.Data;
b.Data.Signal =f_signal;
b.Data.Error  =f_error;
b.Data.Monitor=f_monitor;
for index=1:length(f_axes)
  b.Data.([ 'axis' num2str(index) ]) = f_axes{index};
end

% clear aliases
g = getalias(b); g(1:3) = [];
setalias(b, g);
% make new aliases/axes
setalias(b,'Signal', 'Data.Signal');
setalias(b,'Error',  'Data.Error');
setalias(b,'Monitor','Data.Monitor');

% clear axes
setaxis (b, [], getaxis(b));

for index=1:length(f_axes)
  if index <= length(i_labels)
    b=setalias(b,[ 'axis' num2str(index) ], [ 'Data.axis' num2str(index) ], i_labels{index});
  else
    b=setalias(b,[ 'axis' num2str(index) ], [ 'Data.axis' num2str(index) ]);
  end
  b=setaxis (b, index, [ 'axis' num2str(index) ]);
end
b.Command=cmd;
b = iData_private_history(b, mfilename, a, varargin{:});  
% final check
b = iData(b);

% reset warnings during interp
try
  warning(warn.set);
  warning(warn.get);
catch
  warning(warn);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% private function for interpolation
function f_signal = iData_interp(i_axes, i_signal, f_axes, method)

if isempty(i_signal), f_signal=[]; return; end
if length(i_signal) == numel(i_signal)
  for index=1:length(i_axes)
    x=i_axes{index}; x=x(:); i_axes{index}=x;
  end
end
switch length(i_axes)
case 1    % 1D
  f_signal = interp1(i_axes{1},   i_signal, f_axes{1},   method, 0);
otherwise % nD, n>1
  if length(i_signal) <= 1  % single value ?
    f_signal = i_signal;
    return
  end
  if length(i_signal) == numel(i_signal)  % long vector nD Data set
    if length(i_axes) == 2
      f_signal = griddata(i_axes{:}, i_signal, f_axes{:}, method);
    elseif length(i_axes) == 3
      f_signal = griddata3(i_axes{:}, i_signal, f_axes{:}, method);
    else
      X = cell2mat(i_axes);
      xi= cell2mat(f_axes);
      f_signal = griddatan(X, i_signal, xi, method);
    end
  else
    % f_axes must be an ndgrid result, and monotonic
    f_signal = interpn(i_axes{:}, i_signal, f_axes{:}, method, 0);
  end
end


function b = interp(a, varargin)
% [b...] = interp(s, ...) : interpolate iData object
%
%   @iData/interp function to interpolate data sets
%   This function computes the values of the object 's' interpolated
%   on a new axis grid, which may be specified from an other object, as independent axes,
%   or as a rebinning of the original axes.
%     b=interp(s)    rebin/check 's' on a regular grid.
%     b=interp(s, d) where 'd' is an iData object computes 's' on the 'd' axes.
%     b=interp(s, X1,X2, ... Xn) where 'Xn' are vectors or matrices as obtained 
%                    from ndgrid computes 's' on these axes.
%     b=interp(s, {X1,X2, ... Xn}) is similar to the previous syntax
%     b=interp(s, ..., ntimes) where 'ntimes' is an integer computes new axes for 
%                    interpolation by sub-dividing the original axes ntimes.
%     b=interp(s, ..., 'method') uses specified method for interpolation as one of
%                    linear (default), spline, cubic, or nearest
%     b=interp(s, ..., 'grid') uses meshgrid/ndgrid to determine new axes as arrays
%   Extrapolated data is set to 0 for the Signal, Error and Monitor.
%   For Event data sets, we recommand th use the hist method.
%
% input:  s: object or array (iData)
%         d: single object from which interpolation axes are extracted (iData)
%            or a cell containing axes d={X1,X2, ... Xn}               (cell)
%         X1...Xn: scalar, vectors or matrices specifying axis for 
%            dimensions 1 to ndims(s) (double scalar/vector/matrix)
%         ntimes: original axis sub-division (integer)
% output: b: object or array (iData)
% ex:     a=iData(peaks); b=interp(a, 'grid'); c=interp(a, 2);
%
% Version: $Revision: 1.37 $
% See also iData, interp1, interpn, ndgrid, iData/setaxis, iData/getaxis, iData/hist

% input: option: linear, spline, cubic, nearest
% axes are defined as rank of matrix dimensions
% plot function is plot(y,x,Signal)
% rand(10,20) 10 rows, 20 columns
% pcolor/surf with view(2) shows x=1:20, y=1:10

% handle input iData arrays
if numel(a) > 1
  b = [];
  for index=1:numel(a)
    b = [ b interp(a(index), varargin{:}) ];
  end
  b = reshape(b, size(a));
  return
end

% build new iData object to hold the result
b = copyobj(a);

% object check
if ndims(a) == 0
  iData_private_warning(mfilename,['Object ' inputname(1) ' ' a.Tag ' is empty. Nothing to interpolate.']);
  return
end
% removes warnings during interp
iData_private_warning('enter', mfilename);
warning('off','MATLAB:griddata:DuplicateDataPoints');

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
requires_meshgrid= 0;

% parse varargin to overload defaults and set manually the axes ----------------
axis_arg_index   = 0;
for index=1:length(varargin)
  c = varargin{index};
  if ischar(c) & ~isempty(strfind(c,'grid')) 
    requires_meshgrid=1;
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
  elseif isnumeric(c) & length(c) ~= 1  % set interpolation axes: vector/matrix
    axis_arg_index = axis_arg_index+1;
    if ~isempty(c), f_axes{axis_arg_index} = c; end
  elseif isnumeric(c) & length(c) == 1  % ntimes rebinning (max 10 times)
    ntimes=c;
    if abs(ntimes) > 10, ntimes=10; end
  elseif iscell(c)                      %set interpolation axes: cell(vector/matrix)
    for j1 = 1:length(c(:))
      axis_arg_index = axis_arg_index+1;
      if ~isempty(c{j1}), f_axes{axis_arg_index} = c{j1}; end
    end
  elseif ~isempty(c)
    iData_private_warning(mfilename,['Input argument ' num2str(index) ' of class ' class(c) ' size [' num2str(size(c)) '] is not supported. Ignoring.']);
  end
  clear c
end

b = iData_private_history(b, mfilename, a, varargin{:});
cmd=b.Command;
clear varargin a

% check for method to be valid
if isempty(any(strcmp(method, {'linear','cubic','spline','nearest'})))
  iData_private_warning(mfilename,['Interpolation method ' method ' is not supported. Use: linear, cubic, spline, nearest. Defaulting to linear.']);
  method = 'linear';
end

% check/determine output axes for interpolation --------------------------------

% test axes and decide to call meshgrid if necessary
is_grid=0;
if isvector(b) >= 2 % plot3/event style
  requires_meshgrid=1; 
  if nargin == 1, ntimes=1; end
end

% check for axes for vector: are they vectors, oriented properly ?
if ndims(b) > 1
  for index=1:ndims(b)
    % this axis is a vector, but others are grids: require meshgrid.
    if any(size(f_axes{index}) == 1) & is_grid, requires_meshgrid=1; end 
    try
      % this axis is a grid, others should also be...
      if all(size(f_axes{index}) == size(b)),  is_grid=is_grid+1; end 
    end
  end
end

% trigger regular axis check/rebin when interp(a) called
if nargin == 1 & ~is_grid, ntimes=1; end

% check final axes
for index=1:ndims(b)
  v = f_axes{index}; 
  if isempty(v), v= i_axes{index}; end % no axis specified, use the initial one

  % some of the axes are non consistent or not grid ? re-bin onto vector axes
  if ~is_grid | (requires_meshgrid & is_grid ~= ndims(b)) | ntimes
    % compute the initial axis length
    if isvector(v), a_len = numel(v);
    else            a_len = size( v, index);
    end
    if isvector(b) >= 2 && a_len > prod(size(b))^(1/ndims(b))*2 % event data set
      a_len = prod(size(b))^(1/ndims(b))*2;
    end
    if ntimes > 0 % expand it if requested
      a_len = a_len*ntimes;
    end
    if a_len == 1
      a_len = 2;
    end
    v = v(:); % make all axes vectors
    if numel(v) > numel(i_axes{index})*10 || ntimes > 0 || any(diff(v) == 0)
      v = linspace(min(v), max(v), ceil(a_len) );
    elseif ~issorted(v, 'rows')
      v = unique(v);
    end
  end
  f_axes{index} = v;
  clear v
end

% check if interpolation is indeed required ------------------------------------

% test if interpolation axes have changed w.r.t input object (for possible quick exit)
has_changed = 0;
for index=1:ndims(b)  
  this_i = i_axes{index}; if isvector(this_i), this_i=this_i(:); end
  this_f = f_axes{index}; if isvector(this_f), this_f=this_f(:); end
  if ~isequal(this_i, this_f)
    % length changed ?
    if length(this_i) ~= length(this_f)
      % not same length
      has_changed=1; 
    elseif prod(size(this_i)) ~= prod(size(this_f)) % nb of elements has changed, including matrix axes ?
      has_changed=1; 
    elseif all(abs(this_i(:) - this_f(:)) > 1e-4*abs(this_i(:) + this_f(:))/2)
      % or axis variation bigger than 0.01 percent anywhere
      has_changed=1;
    end
  end
  clear this_i this_f
end

if isvector(b) >=2 % event data set
  b = hist(b, f_axes{:});
  return
end

% get Signal, error and monitor.
i_signal   = get(b,'Signal');

% quick exit is Signal has NaN's (can not interpolate), or axes have not changed
if any(isnan(i_signal(:))), has_changed=1; end
if ~has_changed & (~requires_meshgrid | is_grid), 
  iData_private_warning('exit', mfilename);
  return; 
end

i_class    = class(i_signal); i_signal = double(i_signal);

i_error = getalias(b, 'Error');
if ~isempty(i_error),   
  % check if Error is sqrt(Signal) or a constant
  if strcmp(i_error, 'sqrt(this.Signal)')
    i_error=[];
  elseif isnumeric(i_error) && isscalar(i_error) == 1
    % keep that as  a constant  value
  else
    % else get the value
    i_error  = get(b,'Error');
  end
  i_error    = double(i_error);
end
  
i_monitor = getalias(b, 'Monitor');
if ~isempty(i_monitor),   
  % check if Monitor is 1 or a constant
  if isnumeric(i_monitor) && isscalar(i_monitor) == 1
    % keep that as a constant  value
  else
    % else get the value
    i_monitor  =get(b,'Monitor');
  end
  i_monitor    = double(i_monitor);
end

% do we need to compute a new grid from axis ?
if ndims(b) > 1 && (length(i_signal) ~= numel(i_signal))
  requires_meshgrid = 1;
end

% compute a new grid before interpolating
if requires_meshgrid
  if length(f_axes) == 1 || (ndims(b) == 1 && isvector(b) == 1)
    % nothing to do as we have only one axis, no grid
  else
    % call ndgrid
    [f_axes{1:ndims(b)}] = ndgrid(f_axes{:});  % Warning: may be time consuming
  end
end

% check f_axes vector orientation
for index=1:ndims(b) 
  if isvector(f_axes{index})
    % orient the vector along the dimension
    n = ones(1,ndims(b));
    n(index) = numel(f_axes{index});
    if length(n) == 1, n=[ n 1]; end
    f_axes{index}=reshape(f_axes{index},n);
  end
end

% make sure input axes are monotonic. output axes should be OK.
i_nonmonotonic=0;
for index=1:ndims(b)
  if ~isempty(find(diff(i_axes{index},1,index) < 0))
    i_nonmonotonic=index; break;
  end
end

if i_nonmonotonic
  % transform the initial data into individual points, then interpolate on a regular grid
  for index=1:ndims(b)
    v = i_axes{index}; 
    
    % compute the initial axis length
    if isvector(v), a_len = numel(v);
    else            a_len = size( v, index);
    end
    if isvector(b) >= 2 && a_len > prod(size(b))^(1/ndims(b))*2 % event data set
      a_len = prod(size(b))^(1/ndims(b))*2;
    end
    if a_len == 1
      a_len = 2;
    end
    v = v(:); % make all axes vectors (vectorize initial axes)
    v = linspace(min(v), max(v), ceil(a_len) );
    
    % orient the vector along the dimension
    n = ones(1,ndims(b));
    n(index) = numel(v);
    if length(n) == 1, n=[ n 1]; end
    i_axes_new{index} = reshape(v,n);
    
    % now make the initial axes into event-style data set (all columns)
    v = i_axes{index}; v=v(:); 
    i_axes{index}=v;
    clear v
  end
  i_signal    = iData_interp(i_axes, i_signal(:),  i_axes_new, method);
  if isnumeric(i_error) && length(i_error) > 1, 
    i_error   = iData_interp(i_axes, i_error(:),   i_axes_new, method); 
  end
  if isnumeric(i_monitor) && length(i_monitor) > 1, 
    i_monitor = iData_interp(i_axes, i_monitor(:), i_axes_new, method); 
  end
  i_axes = i_axes_new;
  clear i_axes_new
end

% last test to check if axes have changed
has_changed = 0;
for index=1:ndims(b)    % change to double before interpolation
  i_axes{index}=double(i_axes{index});
  f_axes{index}=double(f_axes{index});
end
for index=1:ndims(b)
  if ~isequal(i_axes{index}, f_axes{index})
    has_changed = 1;
    break
  end
end
if ~has_changed, 
  iData_private_warning('exit', mfilename);
  return; 
end

% interpolation takes place here
f_signal = iData_interp(i_axes, i_signal, f_axes, method);
if isnumeric(i_error) && length(i_error) > 1, 
     f_error = iData_interp(i_axes, i_error,  f_axes, method); 
else f_error = i_error; end
clear i_error
if isnumeric(i_monitor) && length(i_monitor) > 1, 
     f_monitor = iData_interp(i_axes, i_monitor,f_axes, method);
else f_monitor = i_monitor; end
clear i_monitor i_axes

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

% transfer Data and Axes
b.Data.Signal =f_signal;  clear f_signal
b.Data.Error  =f_error;   clear f_error
b.Data.Monitor=f_monitor; clear f_monitor
for index=1:length(f_axes)
  b.Data.([ 'axis' num2str(index) ]) = f_axes{index};
end

% update new aliases, but remove old axes which are numeric (to free memory)
g = getalias(b);
to_remove=[];
for index=4:length(g)
  if any(strcmp(g{index}, b.Alias.Axis))
    if isnumeric(b.Alias.Values{index})
      to_remove=[ to_remove index ];
    end
  end
end
b.Alias.Values(to_remove) = [];
b.Alias.Names(to_remove)  = [];
b.Alias.Labels(to_remove) = [];
setalias(b,'Signal', 'Data.Signal');
setalias(b,'Error',  'Data.Error');
setalias(b,'Monitor','Data.Monitor');

% clear axes
rmaxis (b);

for index=1:length(f_axes)
  if index <= length(i_labels)
    b=setalias(b,[ 'axis' num2str(index) ], [ 'Data.axis' num2str(index) ], i_labels{index});
  else
    b=setalias(b,[ 'axis' num2str(index) ], [ 'Data.axis' num2str(index) ]);
  end
  b=setaxis (b, index, [ 'axis' num2str(index) ]);
end
b.Command=cmd; 
% final check
b = iData(b);

% reset warnings during interp
iData_private_warning('exit', mfilename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% private function for interpolation
function f_signal = iData_interp(i_axes, i_signal, f_axes, method)

if isempty(i_signal), f_signal=[]; return; end
if length(i_signal) == numel(i_signal)
  for index=1:length(i_axes)
    x=i_axes{index}; x=x(:); i_axes{index}=x;
  end
  clear x
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
      f_signal = griddata(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method);
    elseif length(i_axes) == 3
      f_signal = griddata3(i_axes{[2 1 3]}, i_signal, f_axes{[2 1 3]}, method);
    else
      f_signal = griddatan(cell2mat(i_axes), i_signal, cell2mat(f_axes), method);
    end
  else
    % f_axes must be an ndgrid result, and monotonic
    f_signal = interpn(i_axes{:}, i_signal, f_axes{:}, method, 0);
  end
end


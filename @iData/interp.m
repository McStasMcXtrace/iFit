function b = interp(a, varargin)
% [b...] = interp(s, ...) : interpolate iData object
%
%   @iData/interp function to interpolate data sets
%   This function computes the values of the object 's' interpolated
%   on a new axis grid, which may be specified from an other object, as independent axes,
%   or as a rebinning of the original axes.
%     b=interp(s)    rebin/check 's' on a regular grid.
%     b=interp(s, d) where 'd' is an iData object computes 's' on the 'd' axes.
%     b=interp(s, X1,X2, ... Xn) where 'X1...Xn' are vectors or matrices as obtained from ndgrid
%                    computes 's' on these axes.
%     b=interp(s, ..., ntimes) where 'ntimes' is an integer computes new axes for 
%                    interpolation by sub-dividing the original axes ntimes.
%     b=interp(s, ..., 'method') uses specified method for interpolation as one of
%                    linear, spline, cubic, or nearest
%     b=interp(s, ..., 'grid') uses meshgrid/ndgrid to determine new axes as arrays
%   Extrapolated data is set to 0 for the Signal, Error and Monitor.
%
% input:  s: object or array (iData)
%         d: single object from which interpolation axes are extracted (iData)
%            or a cell containing axes d={X1,X2, ... Xn}               (cell)
%         X1...Xn: vectors or matrices specifying axis for dimensions 1 to ndims(s) (double vector/matrix)
%         ntimes: original axis sub-division (integer)
% output: b: object or array (iData)
% ex:     b=interp(a, 'grid');
%
% Version: $Revision: 1.17 $
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
a_axes = cell(1,ndims(a)); a_labels=a_axes;
for index=1:ndims(a)
  [a_axes{index}, a_labels{index}] = getaxis(a, index);  % loads object axes, or 1:end if not defined 
end
for index=ndims(a):length(a.Alias.Axis)
  [dummy, a_labels{index}] = getaxis(a, index);  % additional inactive axes labels (used to create new axes)
end
method='linear';
ntimes=0;

% interpolation axes
i_axes           = a_axes;
axis_arg_index   = 0;
requires_meshgrid= 0;

% parse varargin to overload defaults and set manually the axes
for index=1:length(varargin)
  c = varargin{index};
  if ischar(c) & ~isempty(strfind(c,'grid')) 
    requires_meshgrid=1;
  elseif ischar(c)                      % method (char)
    method = c;
  elseif isnumeric(c) & length(c) > 1   % set interpolation axes: vector/matrix
    axis_arg_index = axis_arg_index+1;
    i_axes{axis_arg_index} = c;
  elseif isnumeric(c) & length(c) == 1  % ntimes rebinning
    ntimes=c;
  elseif iscell(c)
    for j1 = 1:length(c(:))
      axis_arg_index = axis_arg_index+1;
      i_axes{axis_arg_index} = c{j1};
    end
  elseif isa(varargin{index}, 'iData')  % set interpolation axes: get axis from other iData object
    for j1 = 1:ndims(c)
      axis_arg_index = axis_arg_index+1;
      [i_axes{axis_arg_index}, lab] = getaxis(c, j1);
      if ~isempty(lab) & isempty(a_labels{axis_arg_index})
        a_labels{axis_arg_index} = lab;
      end
    end
  else 
    iData_private_warning(mfilename,['Input argument ' num2str(index) ' of class ' class(c) ' size [' num2str(size(c)) '] is not supported. Ignoring.']);
  end
end

% check for method to be valid
if isempty(strmatch(method, {'linear','cubic','spline','nearest'}))
  iData_private_error(mfilename,['Interpolation method ' method ' is not supported. Use: linear, cubic, spline, nearest.']);
end

% test axes and decide to call meshgrid if necessary
is_grid=0;
if isvector(a) >= 2, requires_meshgrid=1; end % plot3/event style
if ndims(a) > 1
  for index=1:ndims(a)
    % this axis is a vector, but others are grids: require meshgrid.
    if any(size(i_axes{index}) == 1)       & is_grid, requires_meshgrid=1; end 
    try
    % this axis is a grid, others should also be...
    if all(size(i_axes{index}) == size(a)),  is_grid=is_grid+1; end 
    catch
    end
  end
end

% trigger regular axis check/rebin when interp(a) called
if nargin == 1 & ~is_grid, ntimes=1; end

if ntimes ~= 0
  % rebin iData object using the smallest axes steps for new axes
  for index=1:ndims(a)
    x = a_axes{index}; x=unique(x); % also makes it a vector
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
      a_len  = min(a_len+1, length(x)*10); % can not reduce or expand more 
                                           % than 10 times each axis
    end
    i_axes{index} = linspace(a_min,a_max,a_len);
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
    v = a_axes{index}; 
    v = unique(v(:));
    i_axes{index} = v;  % vector
  end
end

% test if interpolation axes have changed w.r.t input object
has_changed = 0;
for index=1:ndims(a)  
  if ~isequal(a_axes{index}, i_axes{index})
    % length changed ?
    if length(a_axes{index}) ~= length(i_axes{index})
      % not same length
      has_changed=1; 
    elseif prod(size(a_axes{index})) ~= prod(i_axes{index}) % nb of elements has changed, including matrix axes ?
      has_changed=1; 
    elseif all(abs(a_axes{index} - i_axes{index}) > 1e-4*abs(a_axes{index} + i_axes{index})/2)
      % or axis variation bigger than 0.01 percent anywhere
      has_changed=1;
    end
  end
end

if ~has_changed & (~requires_meshgrid | is_grid), return; end

% get Signal, error and monitor.
a_signal   = get(a,'Signal');
a_class    = class(a_signal); a_signal = double(a_signal);

a_error = getalias(a, 'Error');
if ~isempty(a_error),   
  % check if Error is sqrt(Signal) or a constant
  if strcmp(a_error, 'sqrt(this.Signal)')
    a_error=[];
  elseif isnumeric(a_error) && isscalar(a_error) == 1
    % keep that as  a constant  value
  else
    % else get the value
    a_error  =get(a,'Error');
  end
  a_error    = double(a_error);
end
  
a_monitor = getalias(a, 'Monitor');
if ~isempty(a_monitor),   
  % check if Monitor is 1 or a constant
  if isnumeric(a_monitor) && isscalar(a_monitor) == 1
    % keep that as  a constant  value
  else
    % else get the value
    a_monitor  =get(a,'Monitor');
  end
  a_monitor    = double(a_monitor);
end

% do we need to compute a new grid from axis ?
if ndims(a) > 1 && (length(a_signal) ~= numel(a_signal))
  requires_meshgrid = 1;
end

if requires_meshgrid
  switch ndims(a)
  case 1
    % nothing to do as we have only one axis, no grid
  otherwise
    % calls ndgrid
    % can not use deal as ndgrid uses nargout to set arrays
    toeval='[ ';
    for index=1:ndims(a) 
      toeval=[ toeval 'i_axes{' num2str(index) '}' ];
      if index < ndims(a), toeval=[ toeval ', ' ]; end
    end
    eval([ toeval '] = ndgrid(i_axes{:});' ]);
  end
end

% make sure input axes are monotonic. output axes should be OK.
a_nonmonotonic=0;
for index=1:ndims(a)
  if any(diff(a_axes{index},1,index) < 0)
    a_nonmonotonic=1; break;
  end
end

if a_nonmonotonic
  for index=1:ndims(a)  % apply unique on axes and reorder signal
    a_idx{index}=1:size(a, index);
    [a_axes{index}, a_idx{index}] = unique(a_axes{index});
    if length(a_idx{index}) ~= size(a,index)
      toeval='';
      for j=1:ndims(a), 
        if j ~= index, str_idx{j}=':';
        else           str_idx{j}='a_idx{index}'; end
        if j>1, toeval=[ toeval ',' str_idx{j} ];
        else    toeval=[ str_idx{j} ]; end
      end
      a_signal =eval([ 'a_signal('  toeval ')' ]);
      if isnumeric(a_error) && length(a_error) > 1, 
          try   a_error  =eval([ 'a_error('   toeval ')' ]); 
          catch a_error=[]; end
      end
      if isnumeric(a_error) && length(a_monitor) > 1, 
          try   a_monitor=eval([ 'a_monitor(' toeval ')' ]); 
          catch a_monitor=[]; end
      end
    end
  end
end

i_signal = iData_interp(a_axes, a_signal, i_axes, method);
if isnumeric(a_error) && length(a_error) > 1, 
  i_error  = iData_interp(a_axes, a_error,  i_axes, method); 
else i_error = a_error; end
if isnumeric(a_error) && length(a_monitor) > 1, 
  i_monitor= iData_interp(a_axes, a_monitor,i_axes, method); 
else i_monitor = a_monitor; end

% get back to original Signal class
if ~strcmp(a_class, 'double')
  i_signal = feval(a_class, i_signal);
  i_error  = feval(a_class, i_error);
  i_monitor= feval(a_class, i_monitor);
end

% transfer Data and Axes
Data = a.Data;
Data.Signal =i_signal;
Data.Error  =i_error;
Data.Monitor=i_monitor;
for index=1:length(i_axes)
  Data=setfield(Data,[ 'axis' num2str(index) ], i_axes{index});
end
b.Data = Data;

% clear aliases
g = getalias(b); g(1:3) = [];
setalias(b, g);
% make new aliases/axes
setalias(b,'Signal', 'Data.Signal');
setalias(b,'Error',  'Data.Error');
setalias(b,'Monitor','Data.Monitor');

% clear axes
setaxis (b, [], getaxis(b));
for index=1:length(i_axes)
  b=setalias(b,[ 'axis' num2str(index) ], [ 'Data.axis' num2str(index) ], a_labels{index});
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
function i_signal = iData_interp(a_axes, a_signal, i_axes, method)

if isempty(a_signal), i_signal=[]; return; end
if length(a_signal) == numel(a_signal)
  for index=1:length(a_axes)
    x=a_axes{index}; x=x(:); a_axes{index}=x;
  end
end
switch length(a_axes)
case 1    % 1D
  i_signal = interp1(a_axes{1},   a_signal, i_axes{1},   method, 0);
otherwise % nD, n>1
  if length(a_signal) <= 1  % single value ?
    i_signal = a_signal;
    return
  end
  if length(a_signal) == numel(a_signal)  % long vector nD Data set
    if length(a_axes) == 2
      i_signal = griddata(a_axes{:}, a_signal, i_axes{:}, method);
    elseif length(a_axes) == 3
      i_signal = griddata3(a_axes{:}, a_signal, i_axes{:}, method);
    else
      X = cell2mat(a_axes);
      xi= cell2mat(i_axes);
      i_signal = griddatan(X, a_signal, method);
    end
  else
    % i_axes must be an ndgrid result
    i_signal = interpn(a_axes{:}, a_signal, i_axes{:}, method, 0);
  end
end


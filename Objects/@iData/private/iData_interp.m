function f_signal = iData_interp(i_axes, i_signal, f_axes, method)
% iData_interp: private function for interpolation
% interpolates i_signal(i_axes{}) onto f_axes{} and returns the f_signal
%
% arguments:
%  i_axes:   cell array of initial axes
%  i_signal: initial signal (double array or vector)
%  f_axes:   desired new axes for interpolation (cell array)
%  method:   method used in interpolation
%  f_signal: interpolated new signal (double array)

if isempty(i_signal), f_signal=[]; return; end
if isvector(i_signal) % force axes to be vectors
  for index=1:length(i_axes)
    x=i_axes{index}; x=x(:); i_axes{index}=x;
  end
  clear x
end
if isempty(method), method='linear'; end
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

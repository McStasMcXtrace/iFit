function [f_signal, method] = iData_interp(i_axes, i_signal, f_axes, method)
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

switch length(i_axes)
case 1    % 1D
  X=i_axes{1};
  Y=i_signal;
  [X,I] = unique(X); Y=Y(I);
  f_signal = interp1(X,   Y, f_axes{1},   method, NaN);
otherwise % nD, n>1
  if length(i_signal) <= 1  % single value ?
    f_signal = i_signal;
    return
  end
  if isvector(i_signal)  % long vector nD Data set
    if length(i_axes) == 2
      if ~any(strcmp(method,{'linear','nearest','cubic','v4','natural'})), method='linear'; end
      f_signal = griddata(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method);
      method = 'griddata with vector signal';
    else                       % method: linear or nearest
      if ~any(strcmp(method,{'linear','nearest'})), method='linear'; end
      % i_axes and f_axes must be columns, and cell2mat append them for
      % griddatan
      for index=1:length(i_axes)
        x = i_axes{index}; i_axes{index}=x(:); 
        x = f_axes{index}; f_axes{index}=x(:); clear x;
      end
      f_signal = griddatan(cell2mat(i_axes), i_signal, cell2mat(f_axes), method);
      method = 'griddatan with vector signal';
    end
  else
    % f_axes must be an ndgrid result, and monotonic
    if ~any(strcmp(method,{'linear','nearest','cubic','spline'})), method='linear'; end
    any_vector  = 0;
    any_matrix  = 0;
    for i=1:length(i_axes); 
        v=i_axes{i}; 
        if numel(v) == length(v), i_axes{i}=v(:); any_vector = 1; 
        else any_matrix =1; end; 
    end
    % check if we have mixed axes as vector/matrices -> all go to matrices
    % with repmat. First get the 'matrix' size, which is the one from
    % Signal.
    sz = size(i_signal);
    if any_vector && any_matrix
        for i=1:length(i_axes)
            v=i_axes{i}; 
            if numel(v) == length(v) % isvector
                % orient axis
                this_sz=ones(size(sz)); this_sz(i) = numel(v);
                v=reshape(v, this_sz);
                % replicate it to form a matrix
                this_sz=sz; this_sz(i) = 1;
                v=repmat(v, this_sz);
                i_axes{i} = v;
            end
        end
    end
    % now we can call griddata, else default to interpn
    try
      if length(i_axes) == 2
        if ~any(strcmp(method,{'linear','nearest','cubic','v4','natural'})), method='linear'; end
        f_signal = griddata(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method);
        method = 'griddata with signal';
      else
        if ~any(strcmp(method,{'linear','nearest'})), method='linear'; end
        f_signal = griddatan(i_axes{:}, i_signal, f_axes{:}, method);
        method = 'griddatan with signal';
      end
    catch
      f_signal = interpn(i_axes{:}, i_signal, f_axes{:}, method, NaN);
      method = 'interpn with signal';
    end
  end
end

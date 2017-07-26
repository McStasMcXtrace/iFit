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

f_signal=[];
if isempty(i_signal), 
    return; 
end
if isvector(i_signal) % force axes to be vectors
  for index=1:length(i_axes)
    x=i_axes{index}; x=x(:); i_axes{index}=x;
  end
  clear x
end

% now we test if interpolation is required: are axes the same ?
i_axes2 = i_axes;
if all(cellfun(@(c)numel(c)==length(c), i_axes)) && numel(i_signal) ~= length(i_signal)
  [ i_axes2{:} ] = ndgrid(i_axes{:});
end
same_axes = true;
for index=1:numel(i_axes)
  x1 = i_axes2{index}; x2 = f_axes{index};
  if numel(x1) ~= numel(x2) || any(abs(x1(:) - x2(:)) > 1e-4*x1(:))
    same_axes = false; break
  end
end
if same_axes, 
  f_signal = i_signal; return; 
end

% interpolation is required
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
        if numel(v) == length(v)  % isvector
          sz = ones(size(i_axes));
          sz(i) = numel(v);
          i_axes{i}=reshape(v, sz); 
          any_vector = 1; 
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
    % now we can call interpn, else default to griddata (very slow)
    failed = false;
    f_signal1 = []; 
    try
      f_signal1 = interpn(i_axes{:}, i_signal, f_axes{:}, method, NaN);
      method    = 'interpn with signal';
      f_signal  = f_signal1;
    catch ME
      failed = true;
    end
    % check if we have generated many NaN's
    if ~failed && ~isempty(f_signal1) && length(i_axes) <= 3
      nb_i_nans  = numel(find(isnan(i_signal(:))));
      nb_f_nans1 = numel(find(isnan(f_signal1(:))));
      if nb_f_nans1/numel(f_signal1) > 2*nb_i_nans/numel(i_signal) && nb_f_nans1/numel(f_signal1) > 1e-2
        failed = true;  % twice as many NaN's as before. Try the other interpolant
      end
    end
    
    if failed
      % griddata is much slower than interpn
      if length(i_axes) == 2
        if ~any(strcmp(method,{'linear','nearest','cubic','v4','natural'})), method='linear'; end
        f_signal2 = griddata(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method);
        method = 'griddata with signal';
      elseif length(i_axes) == 3 && exist('griddata3')
        if ~any(strcmp(method,{'linear','nearest'})), method='linear'; end
        f_signal2 = griddata3(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method);
        method = 'griddata3 with signal';
      else
        % try with griddatan, requires to assemble new X,Y and XI
        if numel(i_signal) ~= length(i_signal)
          if all(cellfun(@(c)numel(c)==length(c), i_axes))
            [ i_axes{:} ] = ndgrid(i_axes{:});
          end
          if all(cellfun(@(c)numel(c)==length(c), f_axes))
            [ f_axes{:} ] = ndgrid(f_axes{:});
          end
        end
        % i_axes and f_axes must be columns, and cell2mat append them for
        % griddatan
        sz = size(f_axes{1});
        for index=1:length(i_axes)
          x = i_axes{index}; i_axes{index}=x(:); 
          x = f_axes{index}; f_axes{index}=x(:); clear x;
        end
        f_signal2 = griddatan(cell2mat(i_axes), i_signal, cell2mat(f_axes), method);
        f_signal2 = reshape(f_signal2, sz);
      end
      % get the interpolation which produces less NaN's
      nb_f_nans2 = numel(find(isnan(f_signal2(:))));
      if ~isempty(f_signal1) && ...
        (isempty(f_signal2) || (nb_f_nans1/numel(f_signal1) < nb_f_nans2/numel(f_signal2)))
        f_signal = f_signal1;
      else
        f_signal = f_signal2;
      end
    end
  end
end

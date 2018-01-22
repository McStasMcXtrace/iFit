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
    
  else                   % normal grids/axes
  
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
    
    interpolants = {'iData_interp_interpn','iData_interp_griddata'};
    f_signals ={};
    f_methods ={};
    f_nans    =[];
    % input number of NaN's
    nb_i_nans  = numel(find(isnan(i_signal(:))));
    for index = 1:numel(interpolants)
      try
        [f_signal, f_method] = feval(interpolants{index}, i_axes, i_signal, f_axes, method);
        if ~isempty(f_signal)
          f_signals{end+1}  = f_signal;
          f_methods{end+1}  = f_method;
          f_nans(end+1)     = numel(find(isnan(f_signal(:))));
          if f_nans(end) == 0, break; end % perfect method/no NaN's, we assume it is good
          % less than twice as many NaN's as original data set, or 1%: this is acceptable
          if f_nans(end)/numel(f_signal) < 2*nb_i_nans/numel(i_signal) || f_nans(end)/numel(f_signal) < 1e-2
            break
          end
        end
      catch ME
        disp([ mfilename ': WARNING: failed interpolation with ' interpolants{index} '. Going on.' ]);
      end
    end
    
    f_signal = []; method = []; f_nan = Inf;
    % now we get the best method, which has fewer NaN's
    [~,best] = min(f_nans);
    f_signal = f_signals{best};
    method   = f_methods{best};
    f_nan    = f_nans(best);
  end
end

% ------------------------------------------------------------------------------

function [f_signal, method] = iData_interp_interpn(i_axes, i_signal, f_axes, method)
  % interpn is faster than griddatan, but requires that i_axes are monotonic
  % as obtained from ndgrid.
  switch length(i_axes)
  case 1
    f_signal  = interp1(i_axes{1}, i_signal, f_axes{1}, method, NaN);
  case 2
    f_signal  = interp2(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method, NaN);
  case 3
    f_signal  = interp3(i_axes{[2 1 3]}, i_signal, f_axes{[2 1 3]}, method, NaN);
  otherwise
    f_signal  = interpn(i_axes{:}, i_signal, f_axes{:}, method, NaN);
  end
  method    = [ method ' interpn' ];
  
function [f_signal, method] = iData_interp_griddata(i_axes, i_signal, f_axes, method)
  % griddata is much slower than interpn. 
  % Also, it uses triangulation, which creates fake area when using 
  %   non-meshgrid concave initial data.
  % In order to minimize this effet, we extend the initial data with zeros around
  %
  % NOTE: TriScatteredInterp provides the same results
  
  % pad with 0 around to avoid strange Tri interpolant in concave sets
  [i_axes, i_signal] = iData_interp_safeboxes(i_axes, i_signal);
  
  if length(i_axes) == 2
    if ~any(strcmp(method,{'linear','nearest','cubic','v4','natural'})), method='linear'; end
    f_signal = griddata(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method);
    method = [ method ' griddata'];
  elseif length(i_axes) == 3 && exist('griddata3')
    if ~any(strcmp(method,{'linear','nearest'})), method='linear'; end
    f_signal = griddata3(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method);
    method = [ method ' griddata3'];
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
    f_signal = griddatan(cell2mat(i_axes), i_signal, cell2mat(f_axes), method);
    f_signal = reshape(f_signal, sz);
    method = [ method ' griddatan'];
  end
  
function [i_axes, i_signal] = iData_interp_safeboxes(i_axes, i_signal)
  i_signal = iData_interp_safebox(i_signal);
  for index=1:numel(i_axes)
    i_axes{index} = iData_interp_safebox(i_axes{index});
  end
  
function [x] = iData_interp_safebox(in)
  % surrounds a matrix with a 0's
  x = in;
  s = num2cell(2+size(x));  % this adds 1 every where.
  x(s{:}) = 0;          % fill with NaN's
  if numel(x) == length(x)  % a vector
    x = circshift(x, 1);
  else
    x = circshift(x, ones(1, ndims(x)));
  end

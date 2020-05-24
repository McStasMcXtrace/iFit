function [f_signal, method] = private_interp(i_axes, i_signal, f_axes, method)
% PRIVATE_INTERP private function for interpolation
%   interpolates i_signal(i_axes{}) onto f_axes{} and returns the f_signal
%
% arguments:
%  i_axes:   cell array of initial axes
%  i_signal: initial signal (double array or vector)
%  f_axes:   desired new axes for interpolation (cell array)
%  method:   method used in interpolation
% output:
%  f_signal: interpolated new signal (double array)
%  method:   the interpolation method used

% This file is used by: interp

f_signal=[];
if nargin < 4, method=''; end
if length(i_signal) <= 1
  f_signal = i_signal; return;
end

% interpolation is required
if length(i_axes) == 1    % case: 1D

  if ~any(strcmp(method,{'linear','nearest','cubic','spline','pchip','v5cubic'})), method='linear'; end
  X=i_axes{1};
  Y=i_signal;
  [X,I] = unique(X); Y=Y(I);
  f_signal = interp1(X,   Y, f_axes{1},   method, NaN);

elseif isvector(i_signal)  % case: long vector nD 'event' Data set

  if length(i_axes) == 2    % 2D event (i.e. 2 vector axes)
    if ~any(strcmp(method,{'linear','nearest','cubic','v4'})), method='linear'; end
    i_signal(~isfinite(i_signal)) = 0;
    for index=1:length(i_axes)
        x = i_axes{index}; x(~isfinite(x)) = 0; i_axes{index} = x;
        x = f_axes{index}; x(~isfinite(x)) = 0; f_axes{index} = x;
    end
    warning('off','MATLAB:griddata:DuplicateDataPoints');
    f_signal = griddata(i_axes{[2 1]}, i_signal, f_axes{[2 1]}, method);
    method = 'griddata with vector signal';
  else                       % nD event: method: linear or nearest
    if ~any(strcmp(method,{'linear','nearest'})), method='linear'; end
    % i_axes and f_axes must be columns, and cell2mat append them for
    % griddatan
    for index=1:length(i_axes)
      x = i_axes{index}; i_axes{index}=x(:);
      x = f_axes{index}; f_axes{index}=x(:); clear x;
    end
    warning('off','MATLAB:griddatan:DuplicateDataPoints');
    f_signal = griddatan(cell2mat(i_axes), i_signal, cell2mat(f_axes), method);
    method = 'griddatan with vector signal';
  end

else % normal nD grids/axes

    % f_axes must be an ndgrid result, and monotonic
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

    interpolants = {'interp_interpn','interp_TriScatteredInterp', ...
        'interp_scatteredInterpolant','interp_griddata'};
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
          if (f_nans(end)/numel(f_signal) < 2*nb_i_nans/numel(i_signal) || f_nans(end)/numel(f_signal) < 1e-2)
            break
          end
        end
      catch ME
        warning([ mfilename ': WARNING: failed interpolation with ' interpolants{index} '. Going on.' ]);
        getReport(ME)
      end
    end

    f_signal = []; method = []; f_nan = Inf;
    % now we get the best method, which has fewer NaN's
    [~,best] = min(f_nans);
    if isempty(best) || numel(f_signals) == 1
        best = 1;
    end
    if isempty(f_signals)
        error([ mfilename ': failed all interpolation methods.' ]);
    end
    f_signal = f_signals{best};
    method   = f_methods{best};
    f_nan    = f_nans(best);

end

% ------------------------------------------------------------------------------

function [f_signal, method] = interp_interpn(i_axes, i_signal, f_axes, method)
  % interpn is faster than griddatan, but requires that i_axes are monotonic
  % as obtained from ndgrid.
  %   method: 'nearest', 'linear', or 'spline', 'cubic'.
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

function [f_signal, method] = interp_TriScatteredInterp(i_axes, i_signal, f_axes, method)
  % TriScatteredInterp is recommended in R2010a up to R2013a
  %   method: 'nearest', 'linear', or 'natural'.
  if ~exist('TriScatteredInterp') || length(i_axes) > 3
      f_signal=[]; method=[]; return
  end

  % this interpolant only works with real data
  if ~isreal(i_signal)
    [f_signal, method1] = interp_TriScatteredInterp(i_axes, real(i_signal), f_axes, method);
    f_signal = f_signal +j*interp_TriScatteredInterp(i_axes, imag(i_signal), f_axes, method);
    method = method1;
    return
  end

  % make sure we do not have nan/inf, and use vectors
  i_signal(~isfinite(i_signal)) = 0; i_signal=i_signal(:);
  for index=1:length(i_axes)
      x = i_axes{index}; x(~isfinite(x)) = 0; i_axes{index} = x(:);
      x = f_axes{index}; x(~isfinite(x)) = 0; f_axes{index} = x(:);
  end

  switch length(i_axes)
  case 1
    F = TriScatteredInterp(i_axes{1}, i_signal);
    f_signal = F(f_axes{1});
  case 2
    X = i_axes{2}; Y = i_axes{1};
    Z = i_signal;
    if isvector(X) && isvector(Y) && numel(X)*numel(Y) == numel(Z)
      [X,Y] = meshgrid(X,Y);
      X=X(:); Y=Y(:);
    end
    F = TriScatteredInterp(X,Y,Z);
    XI=f_axes{2}; YI = f_axes{1};
    if isvector(XI) && isvector(YI) && numel(XI) ~= numel(YI)
      [XI,YI] = meshgrid(XI,YI);
      XI=XI(:); YI=YI(:);
    end
    f_signal = F(XI,YI);
  case 3
    F = TriScatteredInterp(i_axes{[2 1 3]}, i_signal);
    f_signal = F(f_axes{[2 1 3]});
  end
  method = 'TriScatteredInterp';


function [f_signal, method] = interp_scatteredInterpolant(i_axes, i_signal, f_axes, method)
  % scatteredInterpolant is recommended above R2013a
  %   method: 'nearest', 'linear', or 'natural'.
  if ~exist('scatteredInterpolant') || ~any(length(i_axes) == [2 3])
      f_signal=[]; method=[]; return
  end

  % make sure we do not have nan/inf, and use vectors
  i_signal(~isfinite(i_signal)) = 0; i_signal=i_signal(:);
  for index=1:length(i_axes)
      x = i_axes{index}; x(~isfinite(x)) = 0; i_axes{index} = x(:);
      x = f_axes{index}; x(~isfinite(x)) = 0; f_axes{index} = x(:);
  end

  switch length(i_axes)
  case 2
    F = scatteredInterpolant(i_axes{[2 1]}, i_signal, method);
    f_signal = F(f_axes{[2 1]});
  case 3
    F = scatteredInterpolant(i_axes{[2 1 3]}, i_signal, method);
    f_signal = F(f_axes{[2 1 3]});
  end
  method = 'scatteredInterpolant';

function [f_signal, method] = interp_griddata(i_axes, i_signal, f_axes, method)
  % griddata is much slower than interpn.
  % Also, it uses triangulation, which creates fake area when using
  %   non-meshgrid concave initial data.
  % In order to minimize this effet, we extend the initial data with zeros around
  %
  % NOTE: TriScatteredInterp provides the same results

  % pad with 0 around to avoid strange Tri interpolant in concave sets
  [i_axes, i_signal] = interp_safeboxes(i_axes, i_signal);

  % make sure we do not have nan/inf
  i_signal(~isfinite(i_signal)) = 0;
  for index=1:length(i_axes)
      x = i_axes{index}; x(~isfinite(x)) = 0; i_axes{index} = x;
      x = f_axes{index}; x(~isfinite(x)) = 0; f_axes{index} = x;
  end

  if length(i_axes) == 2
    if ~any(strcmp(method,{'linear','nearest','cubic','v4','natural'})), method='linear'; end
    % check axes dimensions
    X = i_axes{2};
    Y = i_axes{1};
    Z = i_signal;
    if (isvector(X) && isvector(Y)) || ...
      (all(size(X) == size(Z)) && all(size(Y) == size(Z)))
      warning('off','MATLAB:griddata:DuplicateDataPoints');
      f_signal = griddata(X,Y,Z, f_axes{[2 1]}, method);
      method = [ method ' griddata'];
    else % will not work. return.
      f_signal = [];
      return
    end

  elseif length(i_axes) == 3 && exist('griddata3')
    if ~any(strcmp(method,{'linear','nearest'})), method='linear'; end
    warning('off','MATLAB:griddata3:DuplicateDataPoints');
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
    warning('off','MATLAB:griddatan:DuplicateDataPoints');
    f_signal = griddatan(cell2mat(i_axes), i_signal, cell2mat(f_axes), method);
    f_signal = reshape(f_signal, sz);
    method = [ method ' griddatan'];
  end

% ------------------------------------------------------------------------------
function [i_axes, i_signal] = interp_safeboxes(i_axes, i_signal)
  i_signal = interp_safebox(i_signal);
  for index=1:numel(i_axes)
    i_axes{index} = interp_safebox(i_axes{index});
  end

function [x] = interp_safebox(in)
  % surrounds a matrix with a 0's
  x = in;
  s = num2cell(2+size(x));  % this adds 1 every where.
  x(s{:}) = 0;          % fill with NaN's
  if numel(x) == length(x)  % a vector
    x = circshift(x, 1);
  else
    x = circshift(x, ones(1, ndims(x)));
  end

function s = rotate(a, theta)
% ROTATE Rotate objects around origin.
%   ROTATE(A) rotates around the Signal axis ('vertical') with 180 angular steps.
%   The resulting object dimensionality is increased by 1.
%
%   ROTATE(A, N) rotates using N angular steps in full revolution.
%
%   ROTATE(A, [min max]) and ROTATE(A, [TH0 TH1 TH2 ...]) rotates along the 
%   specified theta range or values. The angles are given in Radians.
%
%
% Example: a=iData(sin(pi*(0:.05:1))); b=rotate(a); plot(b); ...
%   close(gcf); ndims(b) == ndims(a)+1
% Version: $Date$ $Version$ $Author$
% See also iData, iData/camproj

  if nargin < 2,      theta = []; end
  if isempty(theta),  theta = 0; end % will use the default
  if isscalar(theta), 
    if theta<=0, theta = 20; end
    theta = linspace(0, 2*pi, theta);
  elseif length(theta) == 2
    theta = linspace(min(theta), max(theta), 20);
  end

  % handle input iData arrays
  if numel(a) > 1
    s = zeros(iData, numel(a), 1);
    for index=1:numel(a)
      s(index) = feval(mfilename, a(index), theta);
    end
    s = reshape(s, size(a));
    return
  end
  
  % compute vector axes
  s    = copyobj(a);
  s    = interp(s,'grid'); % make sure we use a regular grid (use ndgrid on current object)
  Axes = cell(1, ndims(a)+1);
  for index=1:length(Axes)
    Axes{index} = unique(getaxis(s, index));
  end
  % add Theta
  Axes{end} = theta;
  
  % create a new extended axis set (prepare new object)
  [Axes{:}] = ndgrid(Axes{:});
  
  % now we rotate the X axis, create a new axis (extend dim), and extend Signal
  % with bsxfun
  for index=2:ndims(a)
    s = setaxis(s, index, Axes{index});
  end
  s = setaxis(s, 1,          Axes{1}.*cos(Axes{end}));
  s = setaxis(s, ndims(s)+1, Axes{1}.*sin(Axes{end})); % add axis
  
  % create a '1' vector for the last (new) dimension
  n = ones(1, ndims(a)+1); n(end) = length(theta);
  v = ones(n); % a long vector perpendicular to the current objects
  S = get(s, 'Signal'); if isvector(S) && size(S,1) == 1, S=S'; end
  s = set(s, 'Signal', bsxfun(@times, S, v));

  e = getalias(a, 'Error');
  if ~isnumeric(e), e=get(a, 'Error'); end
  if  isnumeric(e) && length(e) > 1, 
    if isvector(e) && size(e,1) == 1, e=e'; end
    s = set(s, 'Error', bsxfun(@times, e, v));
  end
  clear e
  
  m = getalias(a, 'Monitor');
  if ~isnumeric(m), m=get(a, 'Monitor'); end
  if isnumeric(m) && length(m) > 1, 
    if isvector(m) && size(m,1) == 1, m=m'; end
    s = set(s, 'Monitor', bsxfun(@times, m, v));
  end
  clear m

  s.Command=a.Command;
  s = history(s, mfilename, a, theta);
  

function a = smooth(a, varargin)
% s = smooth(a) : smooth iData objects
%
%   @iData/smooth function to smooth iData objects
%     The smooth method uses a Robust spline smoothing algorithm.
%     You can also use the 'interp' method to smooth objects using e.g. splines,
%     cubic, ... filters.
%
%   Z = SMOOTH(Y) automatically smoothes the uniformly-sampled array Y. Y
%   can be any N-D noisy array (time series, images, 3D data,...). Non
%   finite data (NaN or Inf) are treated as missing values.
%
%   Z = SMOOTH(Y,S) smoothes the array Y using the smoothing parameter S.
%   S must be a real positive scalar. The larger S is, the smoother the
%   output will be. If the smoothing parameter S is omitted (see previous
%   option) or empty (i.e. S = []), it is automatically determined using
%   the generalized cross-validation (GCV) method.
%
%   Z = SMOOTH(Y,W) or Z = SMOOTH(Y,W,S) specifies a weighting array W of
%   real positive values, that must have the same size as Y. Note that a
%   nil weight corresponds to a missing value.
%
%   Robust smoothing
%   ----------------
%   Z = SMOOTHN(...,'robust') carries out a robust smoothing that minimizes
%   the influence of outlying data.
%
%   Reference
%   --------- 
%   Garcia D, Robust smoothing of gridded data in one and higher dimensions
%   with missing values. Computational Statistics & Data Analysis, 2010. 
%   <a href="matlab:web('http://www.biomecardio.com/pageshtm/publi/csda10.pdf')">PDF download</a>
%   http://www.biomecardio.com/matlab/smoothn.html
%
% input:  a: object or array (iData/array of)
% output: s: smoothed data set (iData)
% ex:     c=smooth(a);
%
% Version: $Revision$
% See also iData, iData/interp

% smoothn is in private

% handle input iData arrays
if numel(a) > 1
  for index=1:numel(a)
    a(index) = feval(mfilename, a(index), varargin{:});
  end
  if nargout == 0 & length(inputname(1))
    assignin('caller',inputname(1),a);
  end
  return
end

a = set(a, 'Signal', smoothn(subsref(a,struct('type','.','subs','Signal'), varargin{:})));

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end



function c = convn(a,b)
% CONVN N-dimensional normalised convolution.
%   C = CONVN(A,B) computes the convolution of object A with automatic
%   centering and normalization of the filter B. This is a shortcut for
%   CONV(A,B, 'same pad background center normalize')
%
%   C = CONVN(A, width) convolves A with a normalised Gaussian function
%   which width can be given as a single scalar (same width along all
%   dimensions), or a vector of same length as the object dimension.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/times, estruct/conv, estruct/fft, estruct/xcorr, fconv, fconvn, fxcorr
if nargin ==1
  b=a;
end
if isscalar(b)
  b = [ 1 mean(getaxis(a,1)) double(b) 0]; % use input as a width
  b = gauss(b, getaxis(a,1));
elseif isscalar(a)
  a = [ 1 mean(getaxis(a,1)) double(a) 0]; % use input as a width
  a = gauss(a, getaxis(b,1));
elseif isa(b,'double') && numel(b) == ndims(a)
  p = [];
  g = gauss^(ndims(a));
  ax = {};
  for index=1:ndims(a)
    p = [ p 1 mean(getaxis(a,1)) double(b(index)) 0 ];
    ax{end+1} = getaxis(a, index);
  end
  g = g(p, ax{:});
  b = g;
end
c = binary(a, b, 'conv', 'same pad background center normalize');
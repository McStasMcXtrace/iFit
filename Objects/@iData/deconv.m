function c = deconv(a,b, shape)
% DECONV N-dimensional deconvolution and polynomial division.
%   C = DECONV(A,B) deconvolves object B out of object A.
%   The second argument B should better be centred, normalised, without background.
%
%   C = DECONV(A, width) deconvolves A with a Gaussian function which width
%   can be given as a single scalar (same width along all dimensions),
%   or a vector of same length as the object A dimension.
%
%   C = DECONV(A, B, SHAPE) returns the deconvolution with size and behaviour
%   specified by SHAPE:
%     'full'       returns the full convolution.
%     'same'       (default) returns the central part of the convolution
%                    that is the same size as A.
%     'valid'      returns only those parts of the convolution that are computed
%                    without the zero-padded edges.
%     'pad'        pads the A signal by replicating its starting/ending values
%                    in order to minimize the convolution side effects.
%     'center'     centers the B filter so that convolution does not shift
%                    the A signal.
%     'normalize'  normalizes the B filter so that the convolution does not
%                    change the A signal integral.
%     'background' removes the background from B (subtracts the minimal value)
%     'inverse'    performs an FFT deconvolution (default).
%     'iterative'  performs an iterative deconvolution.
%     'correlation' computes the decorrelation instead of the deconvolution (one
%                    of the FFT's is then conjugated).
%   Default SHAPE is 'deconv same'. Multiple keywords are allowed, for
%   instance 'fft same pad background center normalize'.
%
% Example: a=iData(peaks); b=conv(a); c=deconv(b); ~isempty(c)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/times, iData/convn, iData/fft, iData/xcorr, fconv, fconvn, fxcorr, conv, deconv

if nargin < 3, shape = 'same'; end
if nargin ==1
  b = a;
  shape = 'same pad';
end
if isscalar(b)
  b = [ 1 mean(getaxis(a,1)) double(b) 0]; % use input as a width
  b = gauss(b, getaxis(a,1));
  shape = [ shape ' normalize' ];
elseif isscalar(a)
  a = [ 1 mean(getaxis(b,1)) double(a) 0]; % use input as a width
  a = gauss(a, getaxis(b,1));
  shape = [ shape ' normalize' ];
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
  shape = [ shape ' normalize' ];
end

if ~isempty(strfind(shape,'iter'))
  c = binary(a, b, 'deconv', [ 'iterative ' shape ]);
else
  c = binary(a, b, 'deconv', [ 'deconv ' shape ]);
end

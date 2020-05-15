function c = conv(a,b, shape)
% CONV N-dimensional convolution.
%   C = CONV(A,B) convolves A and B objects. B is often refered as a 'filter'.
%
%   C = CONV(A, width) convolves A with a Gaussian function which width
%   can be given as a single scalar (same width along all dimensions),
%   or a vector of same length as the object dimension.
%
%   C = CONV(A, B, SHAPE) returns the convolution with size and behaviour
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
%     'inverse'    performs an FFT deconvolution.
%     'iterative'  performs an iterative deconvolution.
%     'correlation' computes the correlation instead of the convolution (one
%                    of the FFT's is then conjugated).
%   Default SHAPE is 'same'. Multiple keywords are allowed, for
%   instance 'same pad background center normalize'.
%
%   C = CONV(A) and CONV(A, SHAPE) computes the autoconvolution of A. The default
%   SHAPE is then 'same pad'.
%
%   C = CONV(A, 'tas') convolves with 4D (HKLE) neutron TAS
%   configuration using ResLibCal. Any RESCAL-type parameters are sent to
%   ResLibCal and if a Model exists, it is upgraded with a 4D convolution.
%
% Example: a=estruct(peaks); b=conv(a); ~isempty(b)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/times, estruct/convn, ResLibCal, estruct/fft, estruct/xcorr, fconv, fconvn, fxcorr, conv, deconv

if nargin < 3, shape = 'same'; end
if nargin ==1
  b = a;
  shape = 'same pad';
end

% handle array of objects
if numel(a) > 1
  c = zeros(estruct, numel(a),1);
  for index=1:numel(a)
    c(index) = feval(mfilename, a(index), b, shape);
  end
  c = reshape(c, size(a));
  return
end

% handle input argument types
if isa(b, 'iFunc')
  % we evaluate the Model 'b' with the axes from 'a'
  b = a(b); % call estruct.subsref(iFunc)
elseif ischar(b)
  if strcmp(b, 'tas')
    % convolute the estruct.Model with ResLibCal, overlay parameters from the Data set
    % and provide missing axes from the data set into the Model
    c = ResLibCal(a);
    return
  else
    shape = b; b=a;
  end

elseif isscalar(b)
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

c = binary(a, b, 'conv', shape);

function c = conv(a,b, shape)
% c = conv(a,b) : computes the convolution of iData objects
%
%   @iData/conv function to compute the convolution of data sets (FFT based).
%     A deconvolution mode is also possible.
%     When used with a single scalar value, it is used as a width to build a 
%       gaussian function, with same width along all dimensions
%     When used with a vector of same length as the object dimension, a nD
%       gaussian function with width as vector elements along each diemsions
%
%     The syntax: conv(a, 'tas') upgrades the data object with a neutron TAS
%       configuration using ResLibCal. Any RESCAL-type parameters are sent to 
%       ResLibCal and if a Model exists, it is upgraded with a 4D convolution.
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric or scalar)
%     shape: optional shape of the return value
%          full         Returns the full two-dimensional convolution.
%          same         Returns the central part of the convolution of the same size as a.
%          valid        Returns only those parts of the convolution that are computed
%                       without the zero-padded edges. Using this option, y has size
%                       [ma-mb+1,na-nb+1] when all(size(a) >= size(b)).
%          deconv       Performs an FFT deconvolution.
%          deconv_iter  Performs an iterative deconvolution.
%          pad          Pads the 'a' signal by replicating its starting/ending values
%                       in order to minimize the convolution side effects.
%          center       Centers the 'b' filter so that convolution does not shift
%                       the 'a' signal.
%          normalize    Normalizes the 'b' filter so that the convolution does not
%                       change the 'a' signal integral.
%          background   Remove the background from the filter 'b' (subtracts the minimal value)
%     Default shape is 'same'
%
% output: c: object or array (iData)
% ex:     c=conv(a,b); c=conv(a,b, 'same pad background center normalize');
%
% Version: $Date$
% See also iData, iData/times, iData/convn, iData/fft, iData/xcorr, fconv, fconvn, fxcorr, conv, deconv
if nargin ==1
	b = a;
end
if nargin < 3, shape = 'same'; end

% handle array of objects
if numel(a) > 1
  c = zeros(iData, numel(a),1);
  parfor index=1:numel(a)
    c(index) = feval(mfilename, a(index), b, shape);
  end
  c = reshape(c, size(a));
  return
end

% handle input argument types
if isa(b, 'iFunc')
  % we evaluate the Model 'b' with the axes from 'a'
  b = a(b); % call iData.subsref(iFunc)
elseif strcmp(b, 'tas')
  % convolute the iData.Model with ResLibCal, overlay parameters from the Data set
  % and provide missing axes from the data set into the Model
  c = ResLibCal(a);
  return
  
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

c = iData_private_binary(a, b, 'conv', shape);


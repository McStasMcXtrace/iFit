function c = xcorr(a,b, shape)
% XCORR N-dimensional correlation of objects.
%   C = XCORR(A,B) correlates A and B objects. The correlation is
%   defined as a convolution with the FFT of B conjugated.
%
%   C = XCORR(A, width) correlates A with a Gaussian function which width
%   can be given as a single scalar (same width along all dimensions),
%   or a vector of same length as the object dimension.
%
%   C = CONV(A, B, SHAPE) returns the correlation with size and behaviour
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
%     'inverse'    performs an FFT decorrelation.
%     'iterative'  performs an iterative decorrelation.
%   Default SHAPE is 'same center'. Multiple keywords are allowed, for
%   instance 'same pad background center normalize'.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/times, estruct/convn, estruct/fft, convn, fconv, fconvn
if nargin ==1
  b = a;
end
if nargin < 3, shape = 'same center'; end

c = conv(a, b, [ shape ' correlation' ]);

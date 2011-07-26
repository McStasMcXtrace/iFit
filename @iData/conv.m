function c = conv(a,b, shape)
% c = conv(a,b) : computes the convolution of iData objects
%
%   @iData/conv function to compute the convolution of data sets (FFT based).
%     A deconvolution mode is also possible.
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
%     shape: optional shape of the return value
%          full         Returns the full two-dimensional convolution.
%          same         Returns the central part of the convolution of the same size as a.
%          valid        Returns only those parts of the convolution that are computed
%                       without the zero-padded edges. Using this option, y has size
%                       [ma-mb+1,na-nb+1] when all(size(a) >= size(b)).
%          deconv       Performs an FFT deconvolution.
%          pad          Pads the 'a' signal by replicating its starting/ending values
%                       in order to minimize the convolution side effects
%          center       Centers the 'b' filter so that convolution does not shift
%                       the 'a' signal.
%          normalize    Normalizes the 'b' filter so that the convolution does not
%                       change the 'a' signal integral.
%          background   Remove the background from the filter 'b' (subtracts the minimal value)
%
% output: c: object or array (iData)
% ex:     c=conv(a,b); c=conv(a,b, 'same pad background center normalize');
%
% Version: $Revision: 1.3 $
% See also iData, iData/times, iData/convn, iData/fft, convn, fconv, fconvn
if nargin ==1
	c=conv(a, a); 
	return
end
if nargin < 3, shape = 'same'; end

c = iData_private_binary(a, b, 'conv', shape);


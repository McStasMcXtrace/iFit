function b = ifft(a)
% c = ifft(a,b) : computes the inverse Discrete Fourier transform of iData objects
%
%   @iData/ifft function to compute the inverse Discrete Fourier transform of data sets
%     using the FFT algorithm.
%
% input:  a: object or array (iData)
% output: b: object or array (iData)
% ex:     t=linspace(0,1,1000); 
%         a=iData(t,0.7*sin(2*pi*50*t)+sin(2*pi*120*t)+2*randn(size(t)));
%         c=fft(a); d=ifft(c); plot([ a d ])
%
% Version: $Revision: 1.2 $
% See also iData, iData/fft, conv, convn, CONV2, FILTER, FILTER2, FFT, IFFT

b = fft(a, 'ifft');


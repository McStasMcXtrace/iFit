function b = ifft(a, varargin)
% IFFT Inverse discrete Fourier transform.
%   C = IFFT(A) computes the inverse Discrete Fourier transform of data sets
%   using the FFT algorithm.
%
% Example: t=linspace(0,1,1000); 
%          a=iData(t,0.7*sin(2*pi*50*t)+sin(2*pi*120*t)+2*randn(size(t)));
%          c=fft(a); d=ifft(c); max(abs(a-d(1:1000))) < 1e-10
% Version: $Date$ $Version$ $Author$
% See also iData, iData/fft, iData/conv, FFT, IFFT

b = fft(a, 'ifft', varargin{:});


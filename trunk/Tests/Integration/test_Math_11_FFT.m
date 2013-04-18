function result = Math_11_FFT
  t=linspace(0,1,1000);
  a = iData(t,0.7*sin(2*pi*50*t)+sin(2*pi*120*t)+0.05*randn(size(t)));
  c=fft(a); d=ifft(c);
  if std(abs(a-d)) > 0.3
    result='FAILED';
  else
    result = 'OK  fft ifft';
  end

function result=test_iData_conv

  a=iData(gauss);
  b=convn(a,a);
  c=convn(a,std(a));
  d=xcorr(a,a);
  e=deconv(b,a);
  f=figure;
  h=subplot([a b c d e]);
  if all(0.3 < std([b c d])-std(a)) && all(std([b c d])-std(a) < 0.5)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end

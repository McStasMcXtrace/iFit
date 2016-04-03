function result = test_Models_lorz

  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a,'lorz','','fminimfil');
  b = a(lorz, p);
  fig=figure; plot([ a b ]); close(fig);
  
  if abs(max(abs([ 0.6         1.001      0.0019         0.0068 ])-abs(p))) < 0.01
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end

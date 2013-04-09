function result = Fit_1
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a,'','','fminimfil');
  if abs(max(abs([ 0.61         1.0008      0.0035         0.0001 ])-abs(p))) < 0.01
    result = 'OK  fits(a)';
  else
    result = 'FAILED';
  end 

function result = test_Fit_6_limits
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a, 'gauss', [], 'fminimfil', [ 0.5 0.8 -1 0 ], [ 1 1.2 1 1 ]);
  if abs(max(abs([ 0.62         1.0008      0.0035         0.0001 ])-abs(p))) < 0.01
    result = 'OK  fits(a, ''gauss'', [], ''fminralg'', [ 0.5 0.8 0 0 ], [ 1 1.2 1 1 ]);';
  else
    result = 'FAILED';
  end 

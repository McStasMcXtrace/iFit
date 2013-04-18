function result = Math_12_conv
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  plot([a convn(a) ]);
  a=iData([ ifitpath 'Data/ILL_IN6.dat' ]);
  subplot([a convn(a) ]);
  result = 'OK  conv convn';

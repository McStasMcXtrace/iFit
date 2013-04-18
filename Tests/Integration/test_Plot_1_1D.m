function result = Plot_1_1D
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  plot(a);
  old_mon=getalias(a,'Monitor');
  setalias(a,'Monitor',1);
  figure; plot(a);
  result = 'OK  plot';

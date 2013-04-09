function result = Math_5_proj
  a = iData([ ifitpath 'Data/ILL_IN6.dat' ]);
  xlabel(a, 'Time channel'); % 2nd axis
  ylabel(a, 'Angle channel');% 1st axis
  subplot([ log(a) log(camproj(a)) log(sum(a)) ],'tight');
  subplot([ (a) cumsum(a) ] ,'tight');
  result = 'OK  camproj sum cumsum';

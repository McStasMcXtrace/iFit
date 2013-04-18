function result = test_Plot_9_slices
  a=iData([ ifitpath 'Data/ILL_IN6.dat' ]); 
  plot(a(:,622));                   % extract the object made from channel 622 on second axis, with all columns
  result = 'OK  subsref';

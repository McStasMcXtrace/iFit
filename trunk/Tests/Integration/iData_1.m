function result = iData_1
  a = iData([ ifitpath 'Data/ILL_IN6.dat']);
  get(a);
  a = iData(rand(10));
  if ndims(a) == 2
    result = 'OK  iData([ ifitpath ''Data/ILL_IN6.dat'']);';
  else result = 'FAILED'; end

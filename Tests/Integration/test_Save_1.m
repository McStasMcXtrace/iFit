function result = test_Save_1
  a = iData([ ifitpath 'Data/ILL_IN6.dat']);
  f1= saveas(a,'','pdf');               % save object as a PDF and use object ID as file name
  f2= saveas(a,'MakeItSo','hdf5');      % save object as a HDF5 into specified filename. Extension is appended automatically
  if isempty(f1) || isempty(f2), result='FAILED'; 
  else result = 'OK  saveas';
  end
  try; delete(f1); end
  try; delete(f2); end

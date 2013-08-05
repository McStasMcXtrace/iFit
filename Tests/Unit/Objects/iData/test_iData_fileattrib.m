function result=test_iData_fileattrib

a=iData([ ifitpath 'Data/IRS21360_graphite002_ipg.nxs' ]);


if isstruct(fileattrib(a, 'Signal'))
  result = [ 'OK     ' mfilename ];
else
  result = [ 'FAILED ' mfilename ];
end

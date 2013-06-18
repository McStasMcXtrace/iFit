function result=test_iData_cell

a=iData(peaks);
cell(a);
if iscell(a)
  result = [ 'OK     ' mfilename ];
else
  result = [ 'FAILED ' mfilename ];
end

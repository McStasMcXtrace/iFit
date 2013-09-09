function result = test_iData_smooth

  a = iData([ ifitpath 'Data/Diff_BananaPSD_1314088587.th_y' ]);
  a = a.*(1+0.1*randn(size(a)));
  b = smooth(a);
  c = smooth(a, 'sgolay');
  
  if std(std(double(a-b))) < 0.03
    result = 1;
  else
    result = 0;
  end

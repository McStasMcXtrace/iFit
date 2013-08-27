function result = test_iData_smooth

  a = iData([ ifitpath 'Data/Diff_BananaPSD_1314088587.th_y' ]);
  b = smooth(a);

  result = 1;

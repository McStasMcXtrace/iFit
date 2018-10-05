function result = test_Models_spinwave

  % spin wave
  s=sqw_spinwave('defaults');
  if isempty(s)
    result = [ 'OK     ' mfilename ' SKIPPED: SpinWave not available'];
    return
  end
  try
    S=iData(s, [], 0:.05:4, 0,0, 0:0.5:20);
    result = [ 'OK     ' mfilename ];
  catch
    result = [ 'FAILED ' mfilename ];
  end

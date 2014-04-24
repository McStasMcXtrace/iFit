function result=test_ResLibCal_compute

  out=ResLibCal('default');	% load the default figure/parameters
  out=ResLibCal('compute');
  out=ResLibCal('view2');
  out=ResLibCal('view3');
  out=ResLibCal('geometry');

  ResLibCal('quit');

  if abs(out.resolution.Bragg(4) - 0.5095) < 1e-3
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end


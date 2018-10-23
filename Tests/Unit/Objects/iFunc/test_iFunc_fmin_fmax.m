function result = test_iFunc_fmin_fmax

  model= gauss1;
  fix(model, 'all'); model.Intensity='free';
  model.Intensity=1; model.HalfWidth=.5;
  xlim(model, 'Intensity',[-2 2]);
  p  = fmin(model, [], 'fminpowell');

  p2 = fmax(model, [], 'fminpowell');

  if abs(max(abs([ -2 0 0.5 ])-abs(p))) < 0.1 ...
  && abs(max(abs([  2 0 0.5 ])-abs(p2))) < 0.1
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 

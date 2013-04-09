function result = Math_4_peaks
  a = iData([ ifitpath 'Data/MCA.dat' ]);
  [half_width, center, amplitude, baseline]=peaks(a);
  if length(amplitude) > 100 && length(amplitude) < 120
    result = 'OK  peaks';
  else
    result = 'FAILED';
  end 

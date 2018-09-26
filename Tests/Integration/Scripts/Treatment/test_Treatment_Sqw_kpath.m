function result = test_Treatment_Sqw_kpath

  Sqw   = sqw_cubic_monoatomic('defaults');
  [d,k, fig] = band_structure(Sqw, '','','plot THz');
  close(fig);
  
  if abs(std(d,2)-1.5) < .2
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 
 

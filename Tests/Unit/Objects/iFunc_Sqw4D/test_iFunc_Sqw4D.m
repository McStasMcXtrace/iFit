function result=test_iFunc_Sqw4D

  s=sqw_cubic_monoatomic;
  b=band_structure(s);
  t=thermochemistry(s);
  
  % powder stuff
  q=0.01:0.25:3; w=0:.25:7;
  p=powder(s);
  d=iData(p, [], q, w);
  
  % export to Sqw2D
  file = saveas(d, tempname, 'sqw');
  
  % export to Sqw4D
  d4 = iData(s);
  file4 = saveas(d4, tempname, 'sqw');

  % check value
  if abs(std(d,1) - 0.8) > .1 && exist(file) && exist(file4)
    result = [ 'FAILED ' mfilename ];
  else
    result = [ 'OK     ' mfilename ];
  end
  
  if exist(file),  delete(file);  end
  if exist(file4), delete(file4); end

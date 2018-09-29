function result=test_iFunc_Sqw4D

  s=sqw_cubic_monoatomic;
  b=band_structure(s);
  t=thermochemistry(s);
  
  % powder stuff
  q=0.01:0.25:3; w=0:.25:7;
  p=powder(s);
  d=iData(p, [], q, w);

  % check value
  if abs(std(d,1) - 0.8) > .1
    result = [ 'FAILED ' mfilename ];
  else
    result = [ 'OK     ' mfilename ];
  end

function result=test_iData_vDOS

  s=sqw_cubic_monoatomic;
  gw  = dos(s, 100);
  inc = incoherent(gw, 'T', 300, 'm', 27);
  sinc=deBosify(plus(inc));
  g2 = dos(sinc); g2=g2/trapz(g2);
  
  [~,c1] = std(gw,1);
  [~,c2] = std(g2,1);

  % check value
  if abs(c2 - c1) > .5
    result = [ 'FAILED ' mfilename ];
  else
    result = [ 'OK     ' mfilename ];
  end

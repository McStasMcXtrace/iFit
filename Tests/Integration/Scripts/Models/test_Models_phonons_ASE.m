function result = test_Models_phonons_ASE

  % spin wave
  s=sqw_phonons;
  if isempty(s)
    result = [ 'OK     ' mfilename ' SKIPPED: ASE not available'];
    return
  end
  qh=linspace(0.01,1.5,50);qk=qh; ql=qh'; w=linspace(0.01,50,50);
  f=iData(s,[],qh,qk,ql,w);
  fig=figure; plot3(log(f(:,:,1,:)));
  [w1,c1]=std(f(:,:,1,:),3);
  
  if abs(c1-7.03) < 0.01 && abs(w1-9.8) < 0.1 
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 
  close(fig)

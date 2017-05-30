function result = test_Models_phonons_ASE

  % phonon Al
  s=sqw_phonons('defaults');
  if isempty(s)
    result = [ 'OK     ' mfilename ' SKIPPED: ASE not available'];
    return
  end
  qh=linspace(0.01,1.5,20);qk=qh; ql=qh'; w=linspace(0.01,50,50);
  f=iData(s,[],qh,qk,ql,w);
  fig=figure; plot3(log(f(:,:,1,:)));
  [w1,c1]=std(f(:,:,1,:),3);
  
  if (abs(c1-13) < 3 && abs(w1-9.8) < 0.4) ...
  || (abs(c1-22) < 3 && abs(w1-8) < 0.4)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 
  close(fig)

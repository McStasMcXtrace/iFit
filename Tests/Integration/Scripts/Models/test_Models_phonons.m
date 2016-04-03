function result = test_Models_phonons

  % spin wave
  ac=sqw_sine3d(5);
  qh=linspace(0,1,50);qk=qh; ql=qh'; w=linspace(0.01,10,50);	
  
  
  f=iData(ac,[],qh,qk,ql,w);
  fig=figure; plot3(log(f(:,:,1,:)));
  [w1,c1]=std(f(1,:,:,:));
  
  % perovskite
  s=sqw_vaks('KTaO3'); 
  qh=linspace(0,.5,50);qk=qh; ql=qh; w=linspace(0.01,10,51);
  f=iData(s,[],qh,qk,ql,w); 
  [w2,c2]=std(f(1,:,:,:));
  
  % cubic
  s=sqw_cubic_monoatomic([ 3 3 ]);
  qh=linspace(0,.5,50);qk=qh; ql=qh; w=linspace(0.01,10,51);
  f=iData(s,[],qh,qk,ql,w); 	
  [w3,c3]=std(f(1,:,:,:));
  close(fig);
  
  if c1 < 0.01 && w1 < 0.1 && abs(c2 - 0.15) < .05 && abs(w2 - 0.15) < .05 ...
    && abs(c3 - 0.15) < .05 && abs(w3 - 0.15) < .05
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 

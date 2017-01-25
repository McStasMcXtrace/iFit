function result = test_Models_phonons

  % spin wave
  ac=sqw_sine3d(5);
  qh=linspace(0,1,20);qk=qh; ql=qh'; w=linspace(0.01,10,21);	
  
  
  f=iData(ac,[],qh,qk,ql,w);
  fig=figure; plot3(log(f(:,:,1,:)));
  [w1,c1]=std(f(1,:,:,:));
  
  % perovskite
  s=sqw_vaks('KTaO3'); 
  qh=linspace(0,.5,20);qk=qh; ql=qh'; w=linspace(0.01,10,21);
  f=iData(s,[],qh,qk,ql,w); 
  [w2,c2]=std(f(1,:,:,:));
  
  % cubic
  s=sqw_cubic_monoatomic([ 3 3 ]);
  qh=linspace(0,.5,20);qk=qh; ql=qh'; w=linspace(0.01,10,21);
  f=iData(s,[],qh,qk,ql,w); 	
  [w3,c3]=std(f(1,:,:,:));
  close(fig);

  % acoustopt
  s=sqw_acoustopt(5); 
  qh=linspace(0,.5,50);qk=qh; ql=qh'; w=linspace(0.01,10,50);
  f=iData(s,[],qh,qk,ql,w); 	
  [w4,c4]=std(f(1,:,:,:));
  
  % linquad
  s=sqw_linquad(5); 
  qh=linspace(0,.5,50);qk=qh; ql=qh'; w=linspace(0.01,10,50);
  f=iData(s,[],qh,qk,ql,w); 	
  [w5,c5]=std(f(1,:,:,:));
  
  if c1 < 0.01 && w1 < 0.1 && c2 < 0.01 && w2 < .05 ...
    && c3 < .02 && w3 < .06 && w4 < .15 && w5 < 0.02
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 

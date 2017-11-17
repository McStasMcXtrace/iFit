function f = feval_fast(s)
  % iFunc_Sqw4D: feval_fast: quickly evaluate a 4D Sqw(h=0)
  maxFreq = max(s);
  % evaluate the 4D model onto a mesh filling the Brillouin zone [0:0.5 ]
  qk=linspace(0,0.95,50); qh=0; ql=qk; 
  w =linspace(0.01,maxFreq*1.2,51);
  f =iData(s,[],qh,qk,ql',w);
  f = log(squeeze(f(1,:, :,:)));
  xlabel(f, 'QK [rlu]');
  ylabel(f, 'QL [rlu]');
  zlabel(f, 'Energy [meV]');
  title(f, s.Name);

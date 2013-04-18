function result = Plot_4_overlay_1D
  x=-pi:0.01:pi; a=iData(x,x); 
  a.Error=0;                         % replace default Error=sqrt(Signal) by no-error.
  b=sin(a); c=cos(a); d=exp(-a.*a);  % create new objects by applying operator on the initial linear one
  plot([a b c d]);                   % overlay all objects
  result = 'OK  plot overlay';

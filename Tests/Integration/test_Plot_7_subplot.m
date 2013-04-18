function result = Plot_7_subplot
  x=-pi:0.01:pi; a=iData(x,x); a.Error=0; % replace default Error=sqrt(Signal) by no-error.
  b=sin(a); c=cos(a); d=exp(-a.*a);       % create new objects by applying operator on the initial linear one
  e=iData(flow); f=iData(peaks);          % create 2D and 3D objects
  subplot([a b c d e f]);                 % plot all into a set of separate frames
  result = 'OK  subplot';

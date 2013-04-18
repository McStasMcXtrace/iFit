function result = Math_7_catdog
  x=-pi:0.01:pi; a=iData(x,x); 
  a.Error=0; 
  b=sin(a); c=cos(a); d=exp(-a.*a);
  e=cat(1, [a b c d ]); 
  f=copyobj(e);
  rmaxis(f,1);
  g=cat(2, [a b c d]);
  h=dog(2, g); 
  result = 'OK  cat dog';

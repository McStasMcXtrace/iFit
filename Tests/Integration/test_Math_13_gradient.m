function result = Math_13_gradient
  a=iData(peaks);
  g=gradient(a); 
  if length(g) ~= 2, result='FAILED'; 
  else result = 'OK  gradient';
  end

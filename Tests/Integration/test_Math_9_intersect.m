function result = test_Math_9_intersect
  a = iData(peaks);
  b = copyobj(a);
  a{1} = a{1}+10; a{2} = a{2}+10; 
  a.Signal=a.Signal+5;
  [ai,bi]=intersect(a,b);
  [au,bu]=union(a,b);
  result = 'OK  intersect union';

function result = Math_15_jacobian
  a=iData(peaks); x=linspace(1,2,size(a,1));
  g=jacobian(a, x, [],'half X');
  result = 'OK  jacobian';

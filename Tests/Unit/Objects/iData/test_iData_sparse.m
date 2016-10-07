function result=test_iData_sparse

  a=iData([ ifitpath 'Data/Monitor_GV*']);
  b=hist(a);
  c=b(:,1,:); % sparse only on 1d/2d
  d=sparse(c);
  
  if sum(d,0) == sum(c,0)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end

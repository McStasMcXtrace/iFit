function result = test_Models_mccode

  % test McCode model builder
  y = mccode('defaults');
  % then evaluate the model with raw McCode data
  f = iData(y, [], nan);
  % check dimension and integral
  
  if ndims(f) == 1 && prod(size(f)) == 340 && abs(sum(f,0) - 7.8e4) < 1e3
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 

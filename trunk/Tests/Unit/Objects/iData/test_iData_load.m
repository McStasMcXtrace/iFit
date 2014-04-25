function result=test_iData_load

  tic;
  a = iData([ ifitpath 'Data' ]); % also tests post-loaders scripts
  toc
  
  b=iData([ ifitpath 'Data/30dor.fits' ]);
  c=iData([ ifitpath 'Data/Diff_BananaTheta_1314088587.th' ]);
  
  if all(~isempty([ a ; b ; c ]))
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end

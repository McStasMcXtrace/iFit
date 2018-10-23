function result = test_Models_mccode

  % test McCode model builder
  options.mpirun='none';
  options.dir   = tempname;
  options.ncount = 1e5;
  y = mccode('defaults',options);
  if isempty(y)
    result = [ 'OK     ' mfilename ' (McCode / Instrument not found)' ];
    return
  end
  % set a specific monitor file name
  y.UserData.options.monitor='Diff_BananaTheta*';
  % create 3D view
  [~,fig]=plot(y);
  % then evaluate the model with raw McCode data
  f = iData(y, [], nan);
  % check dimension and integral
  
  if ndims(f) == 1 && prod(size(f)) == 340 && abs(sum(f,0) - 7.8e4) < 1e4
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 
  close(fig);

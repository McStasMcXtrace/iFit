function result = test_Loaders_dir_Data

  tic;
  a = iLoad([ ifitpath 'Data' ]);
  toc
  failed = 0;
  for index=1:length(a)
    if isempty(a{index}),          failed = failed + 1;
    elseif isempty(a{index}.Data), failed = failed + 1;
    end
  end
  if failed
    result = [ 'FAILED ' num2str(failed) '/' num2str(length(a)) ];
  else
    result = [ 'OK     ' mfilename ' (' num2str(numel(a)) ' files)' ];
  end

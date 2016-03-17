function result = test_Loaders_dir_Data

  tic;
  a = iLoad([ ifitpath 'Data' ]);
  toc
  failed = 0;
  for index=1:length(a)
    this = a{index};
    if numel(this) > 1, this=this(1); end
    if isempty(this),          failed = failed + 1;
    elseif isfield(this,'Data') && isempty(this.Data), failed = failed + 1;
    end
  end
  if failed
    result = [ 'FAILED ' num2str(failed) '/' num2str(length(a)) ];
  else
    result = [ 'OK     ' mfilename ' (' num2str(numel(a)) ' files)' ];
  end

function result = Loaders_2
  a = iData([ ifitpath 'Data' ]);
  if length(find(isempty(a))) > 3
    result = [ 'FAILED ' num2str(length(find(isempty(a)))-1) '/' num2str(length(a)) ];
  else
    result = 'OK  load';
  end

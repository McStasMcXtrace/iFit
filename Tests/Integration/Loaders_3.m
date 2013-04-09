function result = Loaders_3
  a = iData(rand(10));
  a = iData(struct('a',1,'b','a string'));
  a = findobj(iData);
  f = figure; peaks;
  a = iData(f);
  result = 'OK  iData(peaks); iData(gcf)';

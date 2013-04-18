function result = Loaders_1
  a = load(iData, [ ifitpath 'Data/ILL_IN6.dat' ]);
  config = iLoad('load config');
  result = 'OK  load; iLoad(''config'')';

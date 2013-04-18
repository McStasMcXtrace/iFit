function result = iData_2_loadarray
  a=load(iData, [ ifitpath 'Data/*.scn']);
  get(a,'Title');
  get(a(1),'Title');
  get(a,'Data.VARIA.A1');
  a(2).Data.VARIA.A1;
  a(3).Data;
  if length(a) ~= 3
    result='FAILED';
  else
    result = 'OK  load(''*.scn'')';
  end

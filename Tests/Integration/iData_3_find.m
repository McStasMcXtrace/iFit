function result = iData_3_find
  a=load(iData, [ ifitpath 'Data/sv1850.scn']);
  [match, field]=findstr(a,'TAS');
  % should return 4 fields
  match{1};
  f=findfield(a,'TAS');
  % should return 2 fields
  if length(match) == 4 && length(f) == 2 
    result = 'OK  findstr; findfield';
  else
    result = 'FAILED';
  end 

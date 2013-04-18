function result = test_iData_4_setalias
  a=load(iData, [ ifitpath 'Data/sv1850.scn']);
  setalias(a,'NewField',42);
  a.NewField = 42;
  set(a,'NewField',42);
  label(a,'NewField','Does god exist ?');
  getalias(a,'NewField');
  setalias(a,'NewField','QH');
  setalias(a,'NewField', '[ this.Data.ZEROS.A1 this.Data.VARIA.A1 ]');
  rmalias(a,'NewField');
  ndims(a); % should be 1
  size(a);  % should be 15 1
  getalias(a,'Signal'); % should be CNTS
  getalias(a,'Error');
  getaxis(a, 1 );
  getaxis(a,'1');
  label(a,1);
  result = 'OK  load;setalias;getaxis;label';

function result = test_Math_6_combine
  a = iData([ ifitpath 'Data/ILL_IN6*.dat' ]);
  a(1)=setalias(a(1),'Monitor', 1);
  a(2)=setalias(a(2),'Monitor', 10);
  b=combine(a);
  c=a(1)+a(2);
  result = 'OK  plus combine';

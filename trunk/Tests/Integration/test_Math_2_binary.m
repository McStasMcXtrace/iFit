function result = test_Math_2_binary
   a = iData([ ifitpath 'Data/ILL_IN6.dat' ]);
   b = iData([ ifitpath 'Data/ILL_IN20.dat' ]); 
   set(a,'Monitor', mean(b.Monitor)/10);
   c=a.*b; c=a./b; c=a.^b;
   c = a+b;
   subplot([ log(a) b log(c) ] ,'tight');
   result = 'OK  plus';

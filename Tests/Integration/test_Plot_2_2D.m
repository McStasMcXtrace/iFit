function result = test_Plot_2_2D
  a=load(iData, [ ifitpath 'Data/ILL_D10.dat' ]);
  plot(a);
  a=iData(peaks);
  plot(a);           % a surface, same as surf(a)
  plot(a,'mesh');    % a wired mesh
  plot(a,'contour'); % contour plot, same as contour(a)
  plot(a,'contourf');% contour plot with filled regions
  plot(a,'surfc');   % a surface with contour plot below
  plot(a,'plot3');   % a surface made of lines side by side
  plot(a,'scatter3');% a surface made of colored points, same as scatter3(a)  
  result = 'OK  surf contour surfc plot3 scatter3';

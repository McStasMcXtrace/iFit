function result = test_Plot_3_3D
  [x,y,z,v]=flow; c=iData(x,y,z,v);
  plot(c);
  plot(c,'surf median');    % plots the c=median(signal) isosurface, same as plot(d) [default]
  plot(c,'surf mean');      % plots the c=mean(signal) isosurface
  plot(c,'surf half');      % plots the c=(max-min)/2 isosurface
  plot(c,'plot3');          % plots a volume rendering with semi-transparent style
  plot(c,'scatter3');       % a set of colored points in space
  slice(c);
  result = 'OK  surf contour surfc plot3 scatter3 slice';

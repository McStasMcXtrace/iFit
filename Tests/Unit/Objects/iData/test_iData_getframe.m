function result=test_iData_getframe

  a=iData(peaks);
  
  gl=opengl;
  opengl software
  if isstruct(getframe(a))
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end
  opengl(gl);

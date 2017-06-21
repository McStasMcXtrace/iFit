function result = test_Process

  pid=Process('ping 127.0.0.1'); 
  r1 = isreal(pid); silent(pid); pause(2); 
  exit(pid); r2=isreal(pid); delete(pid);
  
  if r1 && ~r2
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 
 

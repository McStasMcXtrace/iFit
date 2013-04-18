function result = test_Plot_6_sidebyside
  x=-pi:0.01:pi; a=iData(x,x); 
  a.Error=0;                         % replace default Error=sqrt(Signal) by no-error.
  b=sin(a); c=cos(a); d=exp(-a.*a);  % create new objects by applying operator on the initial linear one
  a{2}=1; b{2}=1.5; c{2}=3; d{2}=5;  % assign a new 2D axis single value to each 1D objects
  plot([a b c d],'surf');            % plot all as a set of lines side by side
  e=cat(2, [a b c d]);               % catenate 1D objects into a 2D object along 2nd axis 
  e{2} = [ 1 1.5 3 5 ];              % asign 2nd axis values in one go
  plot(e,'mesh');                    % plot
  result = 'OK  waterfall cat mesh';

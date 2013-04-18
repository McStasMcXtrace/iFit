function result = test_Plot_5_overlay_2D
  [x,y,z]=peaks; a=iData(x,y*10,z); 
  c=linspace(a,-a+50,10);            % continuously go from 'a' to a '-a+50' in 10 steps
  plot(c);                           % plot all on the same figure
  colormap(c);
  result = 'OK  linspace/waterfall';

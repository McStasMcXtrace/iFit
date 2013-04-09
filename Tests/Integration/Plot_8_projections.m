function result = Plot_8_projections
  a=iData([ ifitpath 'Data/ILL_IN6.dat' ]);                       % import data
  y=a.Data.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF_7; y=y(32:371); a{1}=y; % define the angular axis
  ylabel(a,'Angle [deg]');
  subplot([ log(a) sum(a) trapz(a) camproj(a) ],'axis tight');
  result = 'OK  sum trapz camproj';

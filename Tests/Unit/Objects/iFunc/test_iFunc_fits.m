function result = test_iFunc_fits

  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a,'','','fminimfil');
  
  p2 = fits(a,'',struct(),'fminimfil');
  
  if abs(max(abs([ 0.61 1.0008 0.0035 0.0001 ])-abs(p))) < 0.1 ...
  && isstruct(p2) && abs(max(abs([ 0.61 1.0008 0.0035 0.0001 ])-abs(cell2mat(struct2cell(p2)')))) < 0.1
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end 

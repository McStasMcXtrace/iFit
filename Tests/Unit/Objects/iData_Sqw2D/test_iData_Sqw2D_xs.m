function result=test_iData_Sqw2D_xs

  s=iData_Sqw2D('SQW_coh_lGe.nc');
  sigma = scattering_cross_section(Bosify(symmetrize(s),1235), 14.6);

  % check value
  if abs(sigma - 0.46) > .5
    result = [ 'FAILED ' mfilename ];
  else
    result = [ 'OK     ' mfilename ];
  end

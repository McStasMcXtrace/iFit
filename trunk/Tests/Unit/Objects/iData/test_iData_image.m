function result=test_iData_image

  h=image(iData(peaks),[],[], 'hide axes');
  close(gcf);
  if ~isempty(h)
    result = [ 'OK     ' mfilename ];
  else
    result = [ 'FAILED ' mfilename ];
  end

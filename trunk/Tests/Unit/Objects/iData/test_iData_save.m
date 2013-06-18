function result=test_iData_save

  a=iData(peaks);
  formats={'m','mat','fig','ps','hdf4','jpg','tiff','hdf5','nc','edf','png',...
    'csv','svg','wrl','dat','ply','vtk', 'stl', 'off', 'x3d', 'fits', ...
    'yaml','xml','pdf','eps'};
  for index=1:length(formats)
    f=save(a, 'test', formats{index}); delete(f);
  end
  

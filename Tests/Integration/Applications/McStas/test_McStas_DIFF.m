function result = test_McStas_DIFF

  results = mcstas; % test if installed
  
  if ~isempty(results)
    results = mcstas('templateDIFF', 'RV=1; lambda=2.36; Powder=Na2Ca3Al2F14.laz');
    if ~isa(results, 'iData')
      result = [ 'FAILED ' mfilename  ];
    else
      result = [ 'OK     ' mfilename ];
    end
  else
    result = [ 'OK     ' mfilename ' SKIPPED: McStas/McXtrace not available'];
  end
  
  

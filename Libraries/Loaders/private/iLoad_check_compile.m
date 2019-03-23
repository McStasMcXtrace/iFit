function data = iLoad_check_compile(filename, loader, config)
    if strcmp(loader, 'compile')
      % force compile
      read_anytext('compile');
      read_cbf('compile');
      cif2hkl('compile');
    else
      % make a check of installed MeX/binaries
      try
        config.external.cbf = read_cbf('check');
        disp([ mfilename ': CBF     importer  is functional as ' config.external.cbf ]);
      catch ME
        disp([ mfilename ': CBF     importer  is functional as read_cbf.m (failed mex)' ]);
      end
      try
        config.external.looktxt = read_anytext('check');
        disp([ mfilename ': Text    importer  is functional as ' config.external.looktxt ])
      catch ME
        warning([ mfilename ': Text    importer  is NOT functional' ]);
        disp(getReport(ME))
      end
      try
        config.external.cif2hkl = cif2hkl('check');
        disp([ mfilename ': CIF2HKL converter is functional as ' config.external.cif2hkl ])
      catch ME
        warning([ mfilename ': CIF2HKL converter is NOT functional' ]);
        disp(getReport(ME))
      end
        
    end
    
    if ~isempty(filename)
      data    = iLoad(filename, 'load config', varargin{:});
    else
      data    = config;
    end
    disp([ '% Loaded iLoad format descriptions from ' config.FileName ]);

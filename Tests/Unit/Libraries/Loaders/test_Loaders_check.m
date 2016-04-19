function result = test_Loaders_check

  config = iLoad('check');
  failed = '';
  if ~isfield(config, 'external'), failed = 1;
  else
    if ~isfield(config.external, 'cbf'), failed = [ failed 'cbf ' ]; end
    if ~isfield(config.external, 'looktxt'), failed = [ failed 'looktxt ' ]; end
  % removed FAILURE from cif2hkl as this is not vital
  %  if ~isfield(config.external, 'cif2hkl'), failed = [ failed 'cif2hkl ' ]; end
  end
  
  if ~isempty(failed)
    result = [ 'FAILED ' mfilename ' ' failed ];
  else
    if ~isfield(config.external, 'cif2hkl')
      result = [ 'OK     ' mfilename  ' (no cif2hkl)' ];
    else
      result = [ 'OK     ' mfilename ];
    end
  end

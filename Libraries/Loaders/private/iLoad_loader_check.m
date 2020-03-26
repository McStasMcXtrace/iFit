function data = iLoad_loader_check(file, data, loader)
% ILOAD_LOADER_CHECK make the data pretty looking

  if isempty(data), return; end
  % handle case when a single file generates a data set
  newdata = {};
  if isstruct(data) && length(data)>1
    for index=1:numel(data)
      newdata = { newdata{:} iLoad_loader_check(file, data(index), loader) };
    end
    data = newdata;
    return
  elseif iscellstr(data)
    fprintf(1, 'iLoad: Failed to import file %s with method %s (%s). Got a cell of strings. Ignoring\n', file, loader.name, char(loader.method));
  elseif iscell(data) && numel(data)>1
    for index=1:length(data)
      if ~isempty(data{index})
        this = iLoad_loader_check(file, data{index}, loader);
        newdata = { newdata{:} this };
      end
    end
    data = newdata; % now an array of struct
    return
  elseif iscell(data) && numel(find(~cellfun('isempty', data))) == 1
    data = data{find(~cellfun('isempty', data))};
  end

  name='';
  % check loader
  if isstruct(loader)
    method = loader.method;
    options= loader.options;
    if isfield(loader, 'name'), name=loader.name; end
  else
    method = loader; options=''; 
  end

  if isempty(method), method='iFit/load';
  elseif isa(method, 'function_handle'), method = func2str(method);
  end
  if strcmp(loader, 'variable')
    method='iFit/load';
  end
  
  % check all members of the imported data. Create a 'standardized' structure.
  if isempty(name), name=char(method); end
  if iscell(options), options= cellstr(options{1}); options= [ options{1} ' ...' ]; end
  if ~isfield(data, 'Source')  && ~isfield(data, 'Date') ...
   & ~isfield(data, 'Command') && ~isfield(data,' Data')
    new_data.Data = data;
    % transfer some standard fields as possible
    if isfield(data, 'Source'), new_data.Source= data.Source; end
    if isfield(data, 'Title'),  new_data.Title = data.Title; end
    if isfield(data, 'Date'),   new_data.Date  = data.Date; end
    if isfield(data, 'Label'),  new_data.Label = data.Label; end
    
    data = new_data;
    
  end

  if ~isfield(data, 'Source') && ~isfield(data, 'Filename'),  data.Source = file;
  elseif ( ~isfield(data, 'Source') || ( isfield(data, 'Source') && isempty(data.Source) )) ...
     && isfield(data, 'Filename'), data.Source = data.Filename; end

  if ~isfield(data, 'Title'),   
    [pathname, filename, ext] = fileparts(file);
    if ~strcmp(loader, 'variable'), data.Title  = [ filename ext ' ' name  ];
    else data.Title  = [ filename ext ]; end
  end
  data.Title(data.Title == '%') = '';
  
  if ~isfield(data, 'Date')
    if strcmp(loader, 'variable') data.Date   = clock; 
    else
        try
            d=dir(file); data.Date=d.date; 
        end
    end
  end

  if ~isfield(data, 'Format'),
    if ~strcmp(loader, 'variable'), data.Format  = [ name ' import with Matlab ' char(method) ];  
    else data.Format  = [ 'Matlab ' char(method) ]; end
  end

  if strcmp(loader, 'variable')
    data.Command = [ 'iLoad(' file ', ''' char(method) ''', '''  options '''); % ' method ' method ' ];
  else
    data.Command = [ 'iLoad(''' file ''', ''' char(method) ''', ''' options '''); % ' method ' method' ];
  end

  if ~isfield(data, 'Creator'), data.Creator = [ name ' iFit/load/' char(method) ]; 
  else data.Creator = [ name ' iFit/load/' char(method) ' - ' data.Creator ]; end
  if ~isfield(data, 'User'),
    if isunix
      data.User    = [ getenv('USER') ' running on ' computer ' from ' pwd ];
    else
      data.User    = [ 'User running on ' computer ' from ' pwd ];
    end
  end
  if ~isfield(data, 'Label'), data.Label = ''; end
  if     isfield(loader,'name')   data.Format = loader.name; 
  elseif isfield(loader,'method') data.Format=[ char(loader.method) ' import' ]; 
  elseif ischar(loader)           data.Format=[ loader ' import' ]; end
  data.Loader = loader;
  
end % Load_loader_check

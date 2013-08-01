function attribute = iData_getAttribute(in, fields)
% search if there is a corresponding label (in Headers/Attributes)

  attribute = [];
  
  % 'fields' contain the full path to Signal, e.g. 'Data.<path>'

  % if we use 'Headers' from e.g. read_anytext/looktxt
  % we e.g. relate in.Data.<field path> 
  %       to label in.Data.Headers.<field path>
  if isfield(in.Data, 'Headers')
    if strncmp(fields, 'Data.', length('Data.'))
      fields = [ 'Data.Headers.' fields( (length('Data.')+1):end ) ];
    end
    if isfield(in, fields)
      attribute = get(in, fields); % evaluate to get content of link
    end
  end
  
  % if we use 'Attributes' from e.g. read_hdf/HDF or NetCDF/CDF
  % we e.g. relate in.Data.<group>.<field> 
  %      to label  in.Data.<group>.Attributes.<field>
  if isfield(in.Data, 'Attributes')
    % get group and field names
    lastword_index = find(fields == '.', 2, 'last'); % get the group and the field name
    if isempty(lastword_index)
      lastword = fields; 
      group    = '';
      base     = '';
    elseif isscalar(lastword_index)
      lastword = fields((lastword_index+1):end); 
      group    = fields(1:lastword_index);
      base     = '';
    else 
      lastword = fields( (lastword_index(2)+1):end ); 
      group    = fields( (lastword_index(1)+1):lastword_index(2) ); 
      base     = fields(1:lastword_index(1));
    end
    % we prepend the last word with Attributes. and check for existence
    try
      attribute = get(in, [ base group 'Attributes.' lastword ]); % evaluate to get content of link
    catch
      try
        attribute = get(in, [ base group 'Attributes' ]); % evaluate to get content of link
      end
    end
  end
  
  

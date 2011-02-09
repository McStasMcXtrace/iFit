function filename = iData_private_saveas_hdfnc(a, filename)
% private function to write HDF and NetCDF files
  [fields, types, dims] = findfield(a);
  towrite={};
  for index=1:length(fields(:)) % get all field names
    val=get(a, fields{index});
    if isstruct(val) & length(val) > 1
      val = val(1);
      iData_private_warning(mfilename,[ 'Export of member ' fields{index} ' in object ' inputname(1) ' ' a.Tag ' is a structure array. Selecting only 1st element.' ]);
    if iscellstr(val), 
      val=val(:);
      val(:, 2)={ ';' }; 
      val=val'; 
      val=[ val{1:(end-1)} ];
    end
    if ~isnumeric(val) & ~ischar(val), continue; end
    % make sure field name is valid
    n = fields{index};
    n = n(sort([find(isstrprop(n,'alphanum')) find(n == '_') find(n == '.')]));
    fields{index} = n;
    if strcmp(format,'nc') | strcmp(format,'cdf')
      fields{index} = strrep(fields{index}, '.', '_');
    else
      fields{index} = strrep(fields{index}, '.', filesep);
      if isempty(towrite)
        % initial write wipes out the file
        delete(filename);
        hdf5write(filename, [ filesep 'iData' filesep fields{index} ], val, 'WriteMode', 'overwrite');
      else
        % consecutive calls are appended
        try
          hdf5write(filename, [ filesep 'iData' filesep fields{index} ], val, 'WriteMode', 'append');
        catch
          % object already exists: we skip consecutive write
        end
      end
    end
    if isempty(towrite)
      towrite={ fields{index}, val };
    else
      towrite={ towrite{1:end}, fields{index}, val };
    end
  end
  if strcmp(format,'nc') | strcmp(format,'cdf')
    [path, name, ext] = fileparts(filename);
    cdfwrite(name, towrite);
    filename = [ name '.cdf' ];
  end

end % saveas_hdfnc

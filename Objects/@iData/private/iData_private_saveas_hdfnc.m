function filename = iData_private_saveas_hdfnc(a, filename, format)
% private function to write HDF and NetCDF files
  [fields, types, dims] = findfield(a);
  towrite={};
  for index=1:length(fields(:)) % get all field names
    if isempty(fields{index}), continue; end
    val=get(a, fields{index});
    if strcmp(types{index}, 'hdf5.h5string'), val = char(val.Data); end
    if isstruct(val) & length(val) > 1
      val = val(1);
      iData_private_warning(mfilename,[ 'Export of member ' fields{index} ' in object ' inputname(1) ' ' a.Tag ' is a structure array. Selecting only 1st element.' ]);
    end
    if iscellstr(val), 
      val=val(:);
      val(:, 2)={ ';' }; 
      val=val'; 
      val=[ val{1:(end-1)} ];
    end
    if ~isnumeric(val) && ~ischar(val), continue; end
    if isempty(val),                    continue; end % does not support empty values when writing CDF
    % make sure field name is valid
    n = fields{index};
    n = n(sort([find(isstrprop(n,'alphanum')) find(n == '_') find(n == '.')]));
    fields{index} = n;
    if strcmp(format,'nc') | strcmp(format,'cdf') % NetCDF
      fields{index} = strrep(fields{index}, '.', '_');
      if isempty(towrite)
        towrite={ fields{index}, val };
      else
        towrite={ towrite{1:end}, fields{index}, val };
      end
    else                                          % HDF5
      fields{index} = strrep(fields{index}, '.', filesep);
      if isempty(towrite)
        % initial write wipes out the file
        try
        delete(filename);
        end
        hdf5write(filename, [ filesep 'iData' filesep fields{index} ], val, 'WriteMode', 'overwrite');
        towrite = 'append';
        % write root level attributes: file_name, HFD5_Version, file_time, NeXus_version
      else
        % consecutive calls are appended
        try
          hdf5write(filename, [ filesep 'iData' filesep fields{index} ], val, 'WriteMode', 'append');
          % when object already exists: we skip consecutive write
          % write group attributes: NXclass=NXentry or NXData
          % when encountering Signal, define its attributes: signal=1; axes={axes names}
        end
      end
    end
  end % for
  if strcmp(format,'nc') | strcmp(format,'cdf')
    [p, name, ext] = fileparts(filename);
    filename = fullfile(p, name);
    cdfwrite(filename, towrite); % automatically adds .cdf
    filename = [ filename '.cdf' ];
  end

end % saveas_hdfnc

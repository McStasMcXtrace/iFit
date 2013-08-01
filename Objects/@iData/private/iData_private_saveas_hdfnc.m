function filename = iData_private_saveas_hdfnc(a, filename, format)
% private function to write HDF and NetCDF files

  % export all fields
  [fields, types, dims] = findfield(a);
  
  % this will store the list of fields to write
  towrite={};
  for index=1:length(fields(:)) % scan all field names
    if isempty(fields{index}), continue; end
    
    % get the field value
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
    n = n(isstrprop(n,'alphanum') | n == '_'| n == '.');
    
    % is the field an 'Attribute' ?
    
    % now handle different file formats: HDF5, CDF, NetCDF
    if strcmp(format,'cdf')         % CDF
      % build a list of items to wite in one shot at the end of the routine
      n(n == '.') = '_'; 
      if isempty(towrite)
        towrite={ n, val };
      else
        towrite={ towrite{1:end}, n, val };
      end
      clear val
    elseif strncmpi(format,'hdf',3) % HDF5
      n(n == '.') = filesep; 
      
      if isempty(towrite)
        % initial write wipes out the file
        if ~isempty(dir(filename)), delete(filename); end
        hdf5write(filename, [ filesep 'iData' filesep n ], val, 'WriteMode', 'overwrite');
        towrite = 'append';
        % write root level attributes: file_name, HFD5_Version, file_time, NeXus_version
      else
        % consecutive calls are appended
        try
          hdf5write(filename, [ filesep 'iData' filesep n ], val, 'WriteMode', 'append');
          % when object already exists: we skip consecutive write
          % write group attributes: NXclass=NXentry or NXData
          % when encountering Signal, define its attributes: signal=1; axes={axes names}
        end
      end
    end
  end % for
  
  % write the CDF in one shot
  if strcmp(format,'cdf')
    [p, name, ext] = fileparts(filename);
    filename = fullfile(p, name);
    cdfwrite(filename, towrite); % automatically adds .cdf
    filename = [ filename '.cdf' ];
  end

end % saveas_hdfnc

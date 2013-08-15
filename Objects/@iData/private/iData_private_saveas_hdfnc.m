function filename = iData_private_saveas_hdfnc(a, filename, format)
% private function to write HDF, CDF and NetCDF files
  % format='HDF','CDF','NC' (netCDF)

  % export all fields
  [fields, types, dims] = findfield(a);
  
  if strcmp(format,'cdf')
    [p, name, ext] = fileparts(filename);
    filename       = [ fullfile(p, name) '.cdf' ];
    if exist(filename,'file'), delete(filename); end
  end

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
    n = n(isstrprop(n,'alphanum') | n == '_' | n == '.');
    
    % is the field an 'Attribute' ?
    
    % now handle different file formats: HDF5, CDF, NetCDF

    if strcmpi(format,'cdf')  && (isnumeric(val) || ischar(val))        % CDF
      % we ignore Attributes (if any)
      % n(n == '.') = '_'; 
      % if ~isempty(strfind(n, '.Attributes')), continue; end
      try
        if isempty(towrite) % first access: create file
          cdfwrite(filename, {n, val}, 'WriteMode','overwrite');
          towrite = 'append';
        else
          cdfwrite(filename, {n, val}, 'WriteMode','append');
        end
      catch
        fprintf(1, [mfilename  ': Failed to write ' n ' ' class(val) ' ' mat2str(size(val)) '\n' ]);
      end
    elseif strncmpi(format,'hdf',3) % HDF5
      n(n == '.') = '/'; % group separator 
      
      if isempty(towrite) % first access: create file
        % initial write wipes out the file
        hdf5write(filename, [ '/iData/' n ], val, 'WriteMode', 'overwrite');
        towrite = 'append';
        % write root level attributes: file_name, HFD5_Version, file_time, NeXus_version
      else
        % consecutive calls are appended
        try
          hdf5write(filename, [ '/iData/' n ], val, 'WriteMode', 'append');
          % when object already exists: we skip consecutive write
          % write group attributes: NXclass=NXentry or NXData
          % when encountering Signal, define its attributes: signal=1; axes={axes names}
        end
      end
    elseif strcmpi(format,'nc')
      if isempty(towrite) % first access: create file
        if ~isempty(dir(filename)), delete(filename); end
        ncid = netcdf.create(filename, 'CLOBBER');
        towrite = 'append';
      end
      % create dimensions
      if isvector(val), Dims=length(val); 
      else              Dims=size(val); end
      dimId = [];
      for d=1:length(Dims)
        dimId = [ dimId netcdf.defDim(ncid, [ n '_' num2str(d) ], Dims(d) ) ];
      end
      % get the variable storage class
      c = class(val);
      switch class(val)
        case 'double', t='NC_DOUBLE';
        case 'single', t='NC_FLOAT';
        case 'int8',   t='NC_BYTE';
        case 'char',   t='NC_CHAR';
        case 'int16',  t='NC_SHORT';
        case 'int32',  t='NC_INT';
        % netCDF4 types are converted to NetCDF3
        case 'uint8',  val=int8(val);  t='NC_BYTE';
        case 'uint16', val=int16(val); t='NC_SHORT';
        case 'uint32', val=int32(val); t='NC_INT';
        case 'uint64', val=int32(val); t='NC_INT';
        case 'int64',  val=int32(val); t='NC_INT';
        otherwise, t = ''; continue;
      end
      
      if isempty(t)
        fprintf(1, [mfilename  ': Failed to write ' n ' ' c ' ' mat2str(size(val)) '\n' ]);
        continue
      end
      % create the Variable, and set its value
      varid = netcdf.defVar(ncid, n, t, dimId);
      netcdf.endDef(ncid);
      netcdf.putVar(ncid, varid, val);
      netcdf.reDef(ncid);
    end
  end % for
  
  % close netCDF file
  if strcmpi(format,'nc')
    netcdf.close(ncid);
  end

end % saveas_hdfnc



function data = read_hdf5(filename)
% READ_HDF5 Returns the contents of an HDF5 file as a structure
%   The READ_HDF5 function reads an HDF5 file and returns the contents of
%   that file as the fields of a structure.  Groups are treated as elements
%   of their parent structure.  If the file cannot be opened a -1 is
%   returned.
%
%   data = read_hdf5(filename)
%
% The returned structure has:
% data.<group>.<data>             the main branch holding the HDF4 data sets
% data.<group>.Attributes.<data>  the sub-structure holding all attributes with  
%                 similar structure as the main branch
%
%   Example
%       data=read_hdf5('input.h5');
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_hdf4, read_nc, read_cdf

persistent h5_present

if isempty(h5_present)
  if exist('h5info')
    h5_present = 1;
  else
    h5_present = 0;
  end
end

% read file structure
try
  if h5_present
    data_info = h5info(filename);
  else
    data_info = hdf5info(filename);
    data_info = data_info.GroupHierarchy;
  end
catch
  data = []; % not an HDF file
  return
end

% recursive call
[data, fileID] = getGroup(filename, data_info, h5_present);
% close the hdf5 file after reading all groups
if ~isempty(fileID), H5F.close(fileID); end

% return

end

  % ============================================================
  function [data, fileID] = getGroup(filename, data_info, h5_present, fileID)
  % getGroup: recursively traverse the HDF tree

    if nargin < 4, fileID = []; end
    data = [];
    root = data_info.Name;
    if ~strcmp(root, '/'), root = [ root  '/' ]; end

    % Get group datasets
    nvars   = length(data_info.Datasets);
    for i = 1: nvars
        if h5_present
          val = h5read(filename,[root data_info.Datasets(i).Name]);
          dataID = [root data_info.Datasets(i).Name];
        else
          % hdf5read can stall in R2010a. We use a low-level read
          % val = hdf5read(filename,[data_info.Datasets(i).Name]);
          if isempty(fileID)
            fileID = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
          end
          dataID = [];
          try
            dataID = H5D.open(fileID, [root data_info.Datasets(i).Name]);
          end
          try
            if isempty(dataID), dataID = H5D.open(fileID, [ data_info.Datasets(i).Name]); end
          end
          if ~isempty(dataID)
            val    = H5D.read(dataID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
            H5D.close(dataID);
          else
            disp([ mfilename ': ' filename ': ignoring invalid DataSet ' data_info.Datasets(i).Name ]);
            val = [];
          end
          % H5F.close(fileID); done after reading all groups
        end
        if ~isempty(dataID)
          if strcmp(class(val), 'hdf5.h5string'), val = char(val.Data); end
          if iscellstr(val) && length(val) == 1,  val = char(val); end
          name        = getName(data_info.Datasets(i).Name);
          data.(name) = val; clear val;
        end
          
        % get dataset attributes: group.Attributes.<dataset>.<attribute>
        natts = length(data_info.Datasets(i).Attributes);
        if natts && ~isfield(data,'Attributes'), data.Attributes = []; end
        for j=1:natts
            val = data_info.Datasets(i).Attributes(j).Value;
            if strcmp(class(val), 'hdf5.h5string'), val=char(val.Data); end
            if iscellstr(val) && length(val) == 1,  val = char(val); end
            aname        = getName(data_info.Datasets(i).Attributes(j).Name);
            data.Attributes.(name).(aname) = val;
        end
    end
    
    % handle Links -> Data sets
    if isfield(data_info, 'Links') && ~isempty(data_info.Links)
      nlinks  = length(data_info.Links);
      for i = 1: nlinks
        name = getName(data_info.Links(i).Name);
        if isfield(data_info.Links(i),'Value')
          val  = char(data_info.Links(i).Value);
        elseif isfield(data_info.Links(i),'Target')
          val  = char(data_info.Links(i).Target);
        else
            disp([ mfilename ': ' filename ': ignoring link ' name ]);
            continue; 
        end
        % handle the HDF5 link so that it contains valid names
        val((~isstrprop(val,'alphanum') & val ~= '/') | val == '-' | val == '+') = '_';
        if val(1) == '/', val(1) = ''; end
        val(val == '/') = '.';
        data.(name) = val; % associate the link
        % check if there are associated attributes and link them as well
        [base, group, lastword] = getAttributePath(val);
        % chek if exists: [ base group 'Attributes.' lastword ]
        % or              [ base group 'Attributes' ]
        data.Attributes.(name) = [ base group 'Attributes.' lastword ];
      end
    end

    % get group attributes: group.Attributes.<attribute>
    natts = length(data_info.Attributes);
    if natts && ~isfield(data,'Attributes'), data.Attributes = []; end
    for j=1:natts
        val = data_info.Attributes(j).Value;
        if strcmp(class(val), 'hdf5.h5string'), val=char(val.Data); end
        if iscellstr(val) && length(val) == 1,  val = char(val); end
        name = getName(data_info.Attributes(j).Name);
        data.Attributes.(name) = val;
    end

    % Get each subgroup
    ngroups = length(data_info.Groups);
    for i = 1 : ngroups
      [group, fileID] = getGroup(filename, data_info.Groups(i), h5_present, fileID);
      % assign the name of the group
      name = getName(data_info.Groups(i).Name);
      if ~isempty(name)
        data.(name) = group; 
      end
      clear group;
    end
  end


% ------------------------------------------------------------------------------

function name = getName(location)
% getName: get the HDF5 element Name
  if ~ischar(location) || isempty(location)
    name=[]; return
  end
  [p, name, ext]   = fileparts(location);
  name = [ name ext ];
  name = sanitize_name(name);
end
  
function [base, group, lastword] = getAttributePath(field)
% function to split the entry name into basename, group and dataset
% duplicated from iData_getAttribute

  % get group and field names
  lastword_index = find(field == '.' | field == '/', 2, 'last'); % get the group and the field name
  if isempty(lastword_index)
    lastword = field; 
    group    = '';
    base     = '';                            % Attributes.<field>.
  elseif isscalar(lastword_index)
    lastword = field((lastword_index+1):end); 
    group    = field(1:lastword_index);
    base     = '';                            % <group>.Attributes.<field>
  else 
    lastword = field( (lastword_index(2)+1):end ); 
    group    = field( (lastword_index(1)+1):lastword_index(2) ); 
    base     = field(1:lastword_index(1));    % <basename>.<group>.Attributes.<field>
  end
end

% ------------------------------------------------------------------------------
function name = sanitize_name(name)
  name(~isstrprop(name,'print')) = '';
  name(~isstrprop(name,'alphanum')) = '_';
  if isempty(name), return; end
  if name(1) == '_'
    name = name(find(name ~= '_', 1):end);
  end
  if isstrprop(name(1),'digit')
    name = [ 'x' name ];
  end
end

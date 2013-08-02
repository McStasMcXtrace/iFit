function data = read_hdf5(filename)
%READ_HDF5 Returns the contents of an HDF5 file as a structure
%   The READ_HDF5 function reads an HDF5 file and returns the contents of
%   that file as the fields of a structure.  Groups are treated as elements
%   of their parent structure.  If the file cannot be opened a -1 is
%   returned.
%
% The returned structure has:
% data.<group>.<data>             the main branch holding the HDF4 data sets
% data.<group>.Attributes.<data>  the sub-structure holding all attributes with  
%                 similar structure as the main branch
%
%   Example
%       data=read_hdf5('input.h5');

persistent h5_present

if isempty(h5_present)
  if exist('h5info')
    h5_present = 1;
  else
    h5_present = 0;
  end
end

% read file structure
if h5_present
  data_info = h5info(filename);
else
  data_info = hdf5info(filename);
  data_info = data_info.GroupHierarchy;
end

% recursive call
data = getGroup(filename, data_info);

% return

  % inline function ==============================================================
  function data = getGroup(filename, data_info)
  % getGroup: recursively traverse the HDF tree

    data = [];
    root = data_info.Name;
    if ~strcmp(root, '/'), root = [ root  '/' ]; end

    % Get group datasets
    nvars   = length(data_info.Datasets);
    for i = 1: nvars
        if h5_present
          val = h5read(filename,[root data_info.Datasets(i).Name]);
        else
          val = hdf5read(filename,[data_info.Datasets(i).Name]);
        end
        if strcmp(class(val), 'hdf5.h5string'), val=char(val.Data); end
        name        = getName(data_info.Datasets(i).Name);
        data.(name) = val;
        
        % get dataset attributes: group.Attributes.<dataset>.<attribute>
        natts = length(data_info.Datasets(i).Attributes);
        if natts && ~isfield(data,'Attributes'), data.Attributes = []; end
        for j=1:natts
            val = data_info.Datasets(i).Attributes(j).Value;
            if iscell(data_info.Datasets(i).Attributes(j).Value), val = {1}; end
            if strcmp(class(val), 'hdf5.h5string'), val=char(val.Data); end
            aname        = getName(data_info.Datasets(i).Attributes(j).Name);
            data.Attributes.(name).(aname) = val;
        end
    end

    % get group attributes: group.Attributes.<attribute>
    natts = length(data_info.Attributes);
    if natts && ~isfield(data,'Attributes'), data.Attributes = []; end
    for j=1:natts
        val = data_info.Attributes(j).Value;
        if iscell(data_info.Attributes(j).Value), val = {1}; end
        if strcmp(class(val), 'hdf5.h5string'), val=char(val.Data); end
        name = getName(data_info.Attributes(j).Name);
        data.Attributes.(name) = val;
    end

    % Get each subgroup
    ngroups = length(data_info.Groups);
    for i = 1 : ngroups
      group = getGroup(filename, data_info.Groups(i));
      % assign the name of the group
      name = getName(data_info.Groups(i).Name);
      data.(name) = group;
    end
  end

end

function name = getName(location)
  [p, name, ext]   = fileparts(location);
  name = [ name ext ];
  name(~isstrprop(name,'alphanum') | name == '.' | name == '-' | name == '+') = '_';
  end
  

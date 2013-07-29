function data = read_hdf5(filename)
%READ_HDF5 Returns the contents of an HDF5 file as a structure
%   The READ_HDF5 function reads an HDF5 file and returns the contents of
%   that file as the fields of a structure.  Groups are treated as elements
%   of their parent structure.  If the file file cannot be opened a -1 is
%   returned.
%
%   Example
%       data=read_hdf5('input.h5');

% read file structure
if exist('h5info')
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

data = [];
root = data_info.Name;
if ~strcmp(root, '/'), root = [ root  '/' ]; end

% Get group datasets
nvars   = length(data_info.Datasets);
for i = 1: nvars
    if exist('h5read')
      val = h5read(filename,[root data_info.Datasets(i).Name]);
    else
      val = hdf5read(filename,[data_info.Datasets(i).Name]);
    end
    if strcmp(class(val), 'hdf5.h5string'), val=char(val.Data); end
    [p, name]   = fileparts(data_info.Datasets(i).Name);
    data.(name) = val;
    
    % get dataset attributes
    natts = length(data_info.Datasets(i).Attributes);
    if natts && ~isfield(data,'Attributes'), data.Attributes = []; end
    for j=1:natts
        val = data_info.Datasets(i).Attributes(j).Value;
        if iscell(data_info.Datasets(i).Attributes(j).Value), val = {1}; end
        if strcmp(class(val), 'hdf5.h5string'), val=char(val.Data); end
        [p, aname] = fileparts(data_info.Datasets(i).Attributes(j).Name);
        data.Attributes.(name).(aname) = val;
    end
end

% get group attributes
natts = length(data_info.Attributes);
if natts && ~isfield(data,'Attributes'), data.Attributes = []; end
for j=1:natts
    val = data_info.Attributes(j).Value;
    if iscell(data_info.Attributes(j).Value), val = {1}; end
    if strcmp(class(val), 'hdf5.h5string'), val=char(val.Data); end
    [p, name] = fileparts(data_info.Attributes(j).Name);
    data.Attributes.(name) = val;
end

% Get each subgroup
ngroups = length(data_info.Groups);
for i = 1 : ngroups
  group = getGroup(filename, data_info.Groups(i));
  % assign the name of the group
  [p, name] = fileparts(data_info.Groups(i).Name);
  data.(name) = group;
end



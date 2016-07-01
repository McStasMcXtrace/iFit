function b = iData_private_2lamp(a)
% private function to convert/clean an object to match a LAMP processed 
% workspace
% once created, the LAMP workspace can be saved with:
%  save(b, '', 'hdf5 data');

% LAMP processed workspace:
% root: [HISTORY=string
%        MIN_MAX_VALUES=string (size + min + max)
%        OTHER=comment
%        PARAMETERS=inst parameters (double)
%        SOURCE=instr_name
%        TITLES=string
%        Written_by_LAMP=lamp version
%        creator=OTHER
%        file_name, file_time,user
%    NeXus_version = 4.2.0 ; file_name = ... ; HDF5_Version = 1.6.5 ;
%        file_time = 2010-08-26T09:57:42+01:00]
% entry1: [NX_class = NXentry] group
%   data1: [NX_class = NXdata] group
%     DATA: double array with MONITOR=1.0, Z_VALUE=max signal, long_name=title, signal=1
%     PARAMETERS: string
%     SNAPSHOT: 192x192 icon (uint8)
%     X: double array with axis=1, long_name=xlabel
%     Y: double array with axis=2, long_name=ylabel
%     errors: double array
%

% root level attributes for HDF/NeXus
[majnum minnum relnum] = H5.get_libversion;

b = copyobj(a);

% set the new object empty
b=rmaxis(b);
b=rmalias(b);
b.Data = [];

% search for an existing 'entry' item -------------------------------
all_fields = findfield(a);
match      = all_fields(~cellfun(@isempty, strfind(all_fields, 'entry')));

b.Data.Attributes.NeXus_version= '4.2';
b.Data.Attributes.file_name    = a.Source;
b.Data.Attributes.HDF5_Version = sprintf('%i.%i.%i', majnum, minnum, relnum);
b.Data.Attributes.file_time    = datestr(now);
b.Data.Attributes.HISTORY      = sprintf('%s\n', a.Command{:});
b.Data.Attributes.MIN_MAX_VALUES=[ 'w 1: Float dim=[' num2str(size(a)) '] min=' min(a) ' max=' max(a) ];
b.Data.Attributes.OTHER=char(a);
dummy = findfield(a,'PARAMETERS','exact numeric cache');
if isempty(dummy)
  dummy = findfield(a,'PARAMETERS','cache');
  dummy = get(a,dummy);
  index = find(cellfun(@isstruct, dummy));
  if numel(index) == 1
    dummy = dummy{index};
  else
    dummy = dummy{1};
  end
  b.Data.Attributes.PARAMETERS = class2str(' ', dummy);
else
  b.Data.Attributes.PARAMETERS   = get(a,dummy{end});
end

b.Data.Attributes.SOURCE       = a.User;
b.Data.Attributes.TITLES       = [ a.Title '; ' title(a) ];
b.Data.Attributes.Written_by_LAMP=version(iData);
b.Data.Attributes.creator=b.Data.Attributes.OTHER;

% copy any existing 'entryX'
if ~isempty(match)
  if ~any(strcmp(match, 'Data.entry1'))
    [dummy, index] = min(cellfun('length', match));   % field which path length is smallest
    dummy          = get(a, match{index});            % to avoid sub-structure
    if isstruct(dummy)  % must be a 'group' (structure)
      b.Data.entry1 = dummy; 
      clear dummy
    else
      match=[]; % still request creation of group
    end
  else
    b.Data.entry1 = a.Data.entry1;
  end
end

% create NXentry 'entry1' (if does not exist)
b.Data.entry1.Attributes.NX_class='NXentry';


% create NXdata 'workspace' in 'entry1' ----------------------------
% we overwrite any existing workspace, as only one can exist, and should contain
% the Axes+Signal+Error

% Mantid requires data to be double (64 bits).
b.Data.entry1.data1 = []; % clean previous content if any
b.Data.entry1.data1.Attributes.NX_class = 'NXdata';
b.Data.entry1.data1.errors              = double(get(a, 'Error'));

axes_attr = '';
axis_names = {'X','Y','Z','T','U','V','W'};
for index=1:ndims(a)
  
  lab = getaxis(a, num2str(index));
  if ~ischar(lab), lab=''; end
  
  % Mantid requires that axis_name be 'axisN'
  axis_name = axis_names{index};
  
  if index == 1
    axes_attr = [ axis_name axes_attr ]; % also valid for 1D workspaces
  elseif index == 2
    axes_attr = [ axis_name ',' axes_attr ];
  else
    axes_attr = [ axes_attr ',' axis_name ];
  end
  val = double(getaxis(a, index)); 
  
  b.Data.entry1.data1.(axis_name) = val;
  if ~isempty(label(a, index))
    lab = [ lab ' ' label(a, index) ];
  end
  if ~isempty(lab)
    b.Data.entry1.data1.Attributes.(axis_name).long_name = strtrim(lab); 
  end
  b.Data.entry1.data1.Attributes.(axis_name).axis = index;
end

b.Data.entry1.data1.DATA                  = double(a);

b.Data.entry1.data1.Attributes.DATA.signal= int32(1);
b.Data.entry1.data1.Attributes.DATA.axes  = axes_attr;
b.Data.entry1.data1.Attributes.DATA.MONITORS=1;
if ~isempty(label(a, 0))
  b.Data.entry1.data1.Attributes.DATA.long_name = label(a,0);
  b.Data.entry1.data1.Attributes.DATA.units     = label(a,0);
end

% add the PARAMETERS
match      = all_fields(~cellfun(@isempty, strfind(all_fields, 'PARAMETERS')));
if iscell(dummy)
  % we get the first non Attribute match
  match = match(cellfun(@isempty, strfind(match,'Attributes')));
  [dummy, index] = min(cellfun('length', match));   % field which path length is smallest
  dummy          = get(a, match{index});            % to avoid sub-structure
  b.Data.entry1.data1.PARAMETERS = class2str(' ', dummy);
end

% add the snapshot (getframe)
try
  dummy = getframe(a, [ 192 ]);
  b.Data.entry1.data1.SNAPSHOT = floor(mean(dummy.cdata, 3));
end

